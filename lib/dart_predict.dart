import 'dart:math' as Math;

import 'package:satellite_dart/dopplerFactor.dart';
import 'package:satellite_dart/io.dart';
import 'package:satellite_dart/propagation/gstime.dart';
import 'package:satellite_dart/propagation/propagate.dart';
import 'package:satellite_dart/transforms.dart';
import 'package:satellite_dart/utils/abs.dart';

const xkmper = 6.378137E3; // earth radius (km) wgs84
const astro_unit = 1.49597870691E8; // Astronomical unit - km (IAU 76)
const solar_radius = 6.96000E5; // solar radius - km (IAU 76)
const deg2rad = Math.pi / 180;
const ms2day = 1000 * 60 * 60 * 24; // milliseconds to day
const max_iterations = 250;
const defaultMinElevation = 4; // degrees

class DartPredict {
  _observe(Map<String, dynamic> satrec, List<double>? qth, DateTime start) {
    var eci = this._eci(satrec, start);
    var gmst = this._gmst(start);
    if (!eci.position) {
      return null;
    }
    var geo = eciToGeodetic(eci.position, gmst);

    var solar_vector = this._calculateSolarPosition(start);
    var eclipse = this._satEclipsed(eci.position, solar_vector);

    var track = {
      'eci': eci,
      'gmst': gmst,
      'latitude': geo['latitude'] / deg2rad,
      'longitude': this._boundLongitude(geo['longitude'] / deg2rad),
      'altitude': geo['height'],
      'footprint': 12756.33 * Math.acos(xkmper / (xkmper + geo['height'])),
      'sunlit': !eclipse.eclipsed,
      'eclipseDepth': eclipse.depth / deg2rad
    };

    // If we have a groundstation let's get those additional observe parameters
    if (qth != null && qth.length == 3) {
      var observerGd = {
        'longitude': qth[1] * deg2rad,
        'latitude': qth[0] * deg2rad,
        'height': qth[2]
      };

      var positionEcf = eciToEcf(eci.position, gmst),
          velocityEcf = eciToEcf(eci.velocity, gmst),
          observerEcf = geodeticToEcf(observerGd),
          lookAngles = ecfToLookAngles(observerGd, positionEcf),
          doppler = dopplerFactor(observerEcf, positionEcf, velocityEcf);

      track['azimuth'] = lookAngles['azimuth'] / deg2rad;
      track['elevation'] = lookAngles['elevation'] / deg2rad;
      track['rangeSat'] = lookAngles['rangeSat'];
      track['doppler'] = doppler;
    }

    return track;
  }

  _quickPredict(satrec, qth, start, end) {
    var transit = {};
    var lastel = 0;
    var iterations = 0;

    if (this._badSat(satrec, qth, start)) {
      return null;
    }

    var daynum = this._findAOS(satrec, qth, start);
    if (!daynum) {
      return null;
    }
    transit['start'] = daynum;

    var observed = this._observe(satrec, qth, daynum);
    if (!observed) {
      return null;
    }

    var iel = (observed.elevation as double).round();

    var maxEl = 0, apexAz = 0, minAz = 360, maxAz = 0;

    while (iel >= 0 && iterations < max_iterations && (!end || daynum < end)) {
      lastel = iel;
      daynum = daynum +
          ms2day *
              Math.cos((observed.elevation - 1.0) * deg2rad) *
              Math.sqrt(observed.altitude) /
              25000.0;
      observed = this._observe(satrec, qth, daynum);
      iel = (observed.elevation as double).round();
      if (maxEl < observed.elevation) {
        maxEl = observed.elevation;
        apexAz = observed.azimuth;
      }
      maxAz = Math.max(maxAz, observed.azimuth);
      minAz = Math.min(minAz, observed.azimuth);
      iterations += 1;
    }
    if (lastel != 0) {
      daynum = this._findLOS(satrec, qth, daynum);
    }

    transit['end'] = daynum;
    transit['maxElevation'] = maxEl;
    transit['apexAzimuth'] = apexAz;
    transit['maxAzimuth'] = maxAz;
    transit['minAzimuth'] = minAz;
    transit['duration'] = transit['end'] - transit['start'];

    return transit;
  }

  _badSat(satrec, qth, start) {
    if (qth && !this._aosHappens(satrec, qth)) {
      return true;
    } else if (start && this._decayed(satrec, start)) {
      return true;
    } else {
      return false;
    }
  }

  _aosHappens(satrec, qth) {
    var lin, sma, apogee;
    var meanmo =
        satrec.no * 24 * 60 / (2 * Math.pi); // convert rad/min to rev/day
    if (meanmo == 0) {
      return false;
    } else {
      lin = satrec.inclo / deg2rad;

      if (lin >= 90.0) {
        lin = 180.0 - lin;
      }

      sma = 331.25 * Math.exp(Math.log(1440.0 / meanmo) * (2.0 / 3.0));
      apogee = sma * (1.0 + satrec.ecco) - xkmper;

      if ((Math.acos(xkmper / (apogee + xkmper)) + (lin * deg2rad)) >
          abs(qth[0] * deg2rad)) {
        return true;
      } else {
        return false;
      }
    }
  }

  _decayed(satrec, DateTime start) {
    DateTime satepoch =
        DateTime.utc(satrec.epochyr).add(Duration(days: satrec.epochdays));

    var meanmo =
        satrec.no * 24 * 60 / (2 * Math.pi); // convert rad/min to rev/day
    var drag =
        satrec.ndot * 24 * 60 * 24 * 60 / (2 * Math.pi); // convert rev/day^2

    if (satepoch.millisecondsSinceEpoch +
            ms2day * ((16.666666 - meanmo) / (10.0 * abs(drag))) <
        start.millisecondsSinceEpoch) {
      return true;
    } else {
      return false;
    }
  }

  _findAOS(satrec, qth, start) {
    var current = start;
    var observed = this._observe(satrec, qth, current);
    if (!observed) {
      return null;
    }
    var aostime = 0;
    var iterations = 0;

    if (observed.elevation > 0) {
      return current;
    }
    while (observed.elevation < -1 && iterations < max_iterations) {
      current = current -
          ms2day *
              0.00035 *
              (observed.elevation * ((observed.altitude / 8400.0) + 0.46) -
                  2.0);
      observed = this._observe(satrec, qth, current);
      if (!observed) {
        break;
      }
      iterations += 1;
    }
    iterations = 0;
    while (aostime == 0 && iterations < max_iterations) {
      if (!observed) {
        break;
      }
      if (abs(observed.elevation) < 0.50) {
        // this was 0.03 but switched to 0.50 for performance
        aostime = current;
      } else {
        current = current -
            ms2day *
                observed.elevation *
                Math.sqrt(observed.altitude) /
                530000.0;
        observed = this._observe(satrec, qth, current);
      }
      iterations += 1;
    }
    if (aostime == 0) {
      return null;
    }
    return aostime;
  }

  _findLOS(satrec, qth, start) {
    var current = start;
    var observed = this._observe(satrec, qth, current);
    var lostime = 0;
    var iterations = 0;

    while (lostime == 0 && iterations < max_iterations) {
      if (abs(observed.elevation) < 0.50) {
        // this was 0.03 but switched to 0.50 for performance
        lostime = current;
      } else {
        current = current +
            ms2day *
                observed.elevation *
                Math.sqrt(observed.altitude) /
                502500.0;
        observed = this._observe(satrec, qth, current);
        if (!observed) {
          break;
        }
      }
      iterations += 1;
    }
    return lostime;
  }

  _eci(satrec, DateTime date) {
    return propagate(
      satrec,
      date,
    );
  }

  _gmst(DateTime date) {
    return gstime(
      date.year,
      date.month,
      date.day,
      date.hour,
      date.minute,
      date.second.toDouble(),
    );
  }

  _boundLongitude(longitude) {
    while (longitude < -180) {
      longitude += 360;
    }
    while (longitude > 180) {
      longitude -= 360;
    }
    return longitude;
  }

  _satEclipsed(pos, sol) {
    var sd_earth = Math.asin(xkmper / this._magnitude(pos));
    var rho = this._vecSub(sol, pos);
    var sd_sun = Math.asin(solar_radius / rho.w);
    var earth = this._scalarMultiply(-1, pos);
    var delta = this._angle(sol, earth);

    var eclipseDepth = sd_earth - sd_sun - delta;
    var eclipse;
    if (sd_earth < sd_sun) {
      eclipse = false;
    } else if (eclipseDepth >= 0) {
      eclipse = true;
    } else {
      eclipse = false;
    }
    return {'depth': eclipseDepth, 'eclipsed': eclipse};
  }

  _calculateSolarPosition(start) {
    var time = start / ms2day + 2444238.5; // jul_utc

    var mjd = time - 2415020.0;
    var year = 1900 + mjd / 365.25;
    var T = (mjd + this._deltaET(year) / (ms2day / 1000)) / 36525.0;
    var M = deg2rad *
        ((358.47583 +
                ((35999.04975 * T) % 360) -
                (0.000150 + 0.0000033 * T) * Math.pow(T, 2)) %
            360);
    var L = deg2rad *
        ((279.69668 + ((36000.76892 * T) % 360) + 0.0003025 * Math.pow(T, 2)) %
            360);
    var e = 0.01675104 - (0.0000418 + 0.000000126 * T) * T;
    var C = deg2rad *
        ((1.919460 - (0.004789 + 0.000014 * T) * T) * Math.sin(M) +
            (0.020094 - 0.000100 * T) * Math.sin(2 * M) +
            0.000293 * Math.sin(3 * M));
    var O = deg2rad * ((259.18 - 1934.142 * T) % 360.0);
    var Lsa =
        (L + C - deg2rad * (0.00569 - 0.00479 * Math.sin(O))) % (2 * Math.pi);
    var nu = (M + C) % (2 * Math.pi);
    var R = 1.0000002 * (1 - Math.pow(e, 2)) / (1 + e * Math.cos(nu));
    var eps = deg2rad *
        (23.452294 -
            (0.0130125 + (0.00000164 - 0.000000503 * T) * T) * T +
            0.00256 * Math.cos(O));
    R = astro_unit * R;

    return {
      'x': R * Math.cos(Lsa),
      'y': R * Math.sin(Lsa) * Math.cos(eps),
      'z': R * Math.sin(Lsa) * Math.sin(eps),
      'w': R
    };
  }

  _deltaET(year) {
    return 26.465 +
        0.747622 * (year - 1950) +
        1.886913 * Math.sin((2 * Math.pi) * (year - 1975) / 33);
  }

  _vecSub(v1, v2) {
    var vec = {'x': v1.x - v2.x, 'y': v1.y - v2.y, 'z': v1.z - v2.z};
    vec['w'] = this._magnitude(vec);
    return vec;
  }

  _scalarMultiply(k, v) {
    return {
      'x': k * v.x,
      'y': k * v.y,
      'z': k * v.z,
      'w': v.w ? abs(k) * v.w : null
    };
  }

  _magnitude(v) {
    return Math.sqrt(Math.pow(v.x, 2) + Math.pow(v.y, 2) + Math.pow(v.z, 2));
  }

  _angle(v1, v2) {
    var dot = (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
    return Math.acos(dot / (this._magnitude(v1) * this._magnitude(v2)));
  }

  observe(tle, qth, start) {
    var tles = tle.split('\n');
    var satrec = twoline2satrec(tles[1], tles[2]);

    if (this._badSat(satrec, qth, start)) {
      return null;
    }

    return this._observe(satrec, qth, start);
  }

  observes(tle, qth, DateTime start, DateTime end, interval) {
    var tles = tle.split('\n');
    var satrec = twoline2satrec(tles[1], tles[2]);

    if (this._badSat(satrec, qth, start)) {
      return null;
    }

    var observes = [], observed;
    var iterations = 0;
    while (start.millisecondsSinceEpoch < end.millisecondsSinceEpoch &&
        iterations < max_iterations) {
      observed = this._observe(satrec, qth, start);
      if (!observed) {
        break;
      }
      observes.add(observed);
      start.add(interval);
      iterations += 1;
    }

    return observes;
  }

  transits(tle, qth, DateTime start, DateTime end, minElevation, maxTransits) {
    if (!minElevation) {
      minElevation = defaultMinElevation;
    }

    if (!maxTransits) {
      maxTransits = max_iterations;
    }

    var tles = tle.split('\n');
    var satrec = twoline2satrec(tles[1], tles[2]);
    if (this._badSat(satrec, qth, start)) {
      return [];
    }

    var time = start.millisecondsSinceEpoch;
    var transits = [];
    var iterations = 0;

    while (iterations < max_iterations && transits.length < maxTransits) {
      var transit = this._quickPredict(satrec, qth, start, end);
      if (!transit) {
        break;
      }
      if (transit.end > end.millisecondsSinceEpoch) {
        break;
      }
      if (transit.end > start.millisecondsSinceEpoch &&
          transit.maxElevation > minElevation) {
        transits.add(transit);
      }
      time = transit.end + 60 * 1000;
      iterations += 1;
    }

    return transits;
  }

  transitSegment(tle, qth, DateTime start, DateTime end) {
    var tles = tle.split('\n');
    var satrec = twoline2satrec(tles[1], tles[2]);
    if (this._badSat(satrec, qth, start)) {
      return [];
    }

    return this._quickPredict(satrec, qth, start, end);
  }
}

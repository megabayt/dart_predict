# dart-predict
[![Pub Version](https://img.shields.io/pub/v/dart_predict?logo=dart&logoColor=white)](https://pub.dev/packages/dart_predict/)
[![Dart SDK Version](https://badgen.net/pub/sdk-version/dart_predict)](https://pub.dev/packages/dart_predict/)
[![License](https://img.shields.io/github/license/megabayt/dart_predict)](https://raw.githubusercontent.com/megabayt/dart_predict/master/LICENSE)
[![Pub popularity](https://badgen.net/pub/popularity/dart_predict)](https://pub.dev/packages/dart_predict/score)
[![GitHub popularity](https://img.shields.io/github/stars/megabayt/dart_predict?logo=github&logoColor=white)](https://github.com/megabayt/dart_predict/stargazers)

A simple prediction module that implements:

Simple linear regression for basic forecasting http://en.wikipedia.org/wiki/Simple_linear_regression

Moving average to predict next value and remove noize http://en.wikipedia.org/wiki/Moving_average

## Installation

	npm install predict

## Usage
	
	var predict = require('predict');

	//provide input to the linear regression the first array is y and the second is x
	var lr = predict.linearRegression([0,2,4,6],[0,1,2,3]);
	//predict the future result for the value x=4 
	var result = lr.predict(4);	
	//result will be y=8

	//perform moving average with a buffer of size 10 by default 
	var ma = predict.movingAverage();
    ma.pushValues([2,4,6,8,10]);
    //predict next value based on average 2+4+6+8+10 = 30 and 30/5 = 6
    var result = ma.predictNextValue();
    //result will be 6

	//perform moving average with a buffer of size 2 
	var ma = predict.movingAverage(2);
    ma.pushValues([2,4,6,8,10]);
    //predict next value based on average 8+10 = 18 and 18/2 = 9
    var result = ma.predictNextValue();
    //result will be 9
    //as you add values, the average changes
    ma.pushValues([20]);
    result = ma.predictNextValue();
    //will be 15

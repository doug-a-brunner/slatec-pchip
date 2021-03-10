# slatec-pchip

Shape-preserving piecewise cubic Hermite interpolation in pure JS, ported from SLATEC

## Usage

```js
const { Pchip, pchip } = require('slatec-pchip');

const x = [1, 2, 3, 4, 5];
const y = [1, 7, 11, 14, 28];

const interp = new Pchip(x, y);
interp.evaluate(3); // 11
interp.evaluate(4.2); // 15.464470588

pchip(x, y, 3); // 11
pchip(x, y, [1, 4.2]) // [1, 15.464470588]
```

## Details

This is a library for cubic interpolation written in pure JS. It should work as a drop-in replacement for [`interpo`](https://www.npmjs.com/package/interpo), and also allows interpolating values on-the-fly without creating a `Pchip` instance.


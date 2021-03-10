/**
 * Port of SLATEC function DPCHST
 * @param arg1
 * @param arg2
 * @returns 1 if arg1 is the same sign as arg2, 0 if either of them is zero, -1 if arg1 is a different sign from arg2
 */
export function dpchst(arg1: number, arg2: number) {
    if (arg1 === 0 || arg2 === 0) {
        return 0;
    }
    if (arg1 > 0) {
        return arg2 > 0 ? 1 : -1;
    }
    return arg2 > 0 ? -1 : 1;
}

/**
 * Port of SLATEC function DPCHIM: calculates dy/dx at each interpolating point to produce a shape-preserving spline
 * @param x Independent variable values
 * @param y Dependent variable values
 * @returns Derivative dy/dx
 */
export function dpchim(x: number[], y: number[]) {
    if (x.length !== y.length) {
        throw new Error('input array lengths must match');
    }
    const n = x.length;
    if (n < 2) {
        throw new Error('number of data points less than two');
    }
    for (let i = 1; i < n; ++i) {
        if (x[i] <= x[i - 1]) {
            throw new Error('x-array not strictly increasing');
        }
    }
    if (n === 2) {
        const deriv = (y[1] - y[0]) / (x[1] - x[0]);
        return [deriv, deriv];
    }
    const d = new Array(n);
    let h1 = x[1] - x[0];
    let del1 = (y[1] - y[0]) / h1;
    let h2 = x[2] - x[1];
    let del2 = (y[2] - y[1]) / h2;

    // set d[0] via non-centered three-point formula, adjusted to be shape-preserving
    let hsum = h1 + h2;
    let w1 = (h1 + hsum) / hsum;
    let w2 = -h1 / hsum;
    d[0] = w1 * del1 + w2 * del2;
    if (dpchst(d[0], del1) < 0) {
        d[0] = 0;
    }
    else if (dpchst(del1, del2) < 0) {
        // need do this check only if monotonicity switches
        const dmax = 3 * del1;
        if (Math.abs(d[0]) > Math.abs(dmax)) {
            d[0] = dmax;
        }
    }

    // loop through interior points
    for (let i = 1; i < n - 1; ++i) {
        if (i > 1) {
            h1 = h2;
            h2 = x[i + 1] - x[i];
            hsum = h1 + h2;
            del1 = del2;
            del2 = (y[i + 1] - y[i]) / h2;
        }
        d[i] = 0;
        if (dpchst(del1, del2) > 0) {
            // use Brodlie modification of Butland formula
            const hsumt3 = hsum * 3;
            w1 = (hsum + h1) / hsumt3;
            w2 = (hsum + h2) / hsumt3;
            const dmax = Math.max(Math.abs(del1), Math.abs(del2));
            const dmin = Math.min(Math.abs(del1), Math.abs(del2));
            const drat1 = del1 / dmax;
            const drat2 = del2 / dmax;
            d[i] = dmin / (w1 * drat1 + w2 * drat2);
        }
        else {
            d[i] = 0; // set d[i] = 0 unless data are strictly monotonic
        }
    }

    // set d[n - 1] via non-centered three-point formula, adjusted to be shape-preserving
    w1 = -h2 / hsum;
    w2 = (h2 + hsum) / hsum;
    d[n - 1] = w1 * del1 + w2 * del2;
    if (dpchst(d[n - 1], del2) < 0) {
        d[n - 1] = 0;
    }
    else if (dpchst(del1, del2) < 0) {
        // need do this check only if monotonicity switches
        const dmax = 3 * del2;
        if (Math.abs(d[n - 1]) > Math.abs(dmax)) {
            d[n - 1] = dmax;
        }
    }

    return d;
}

/**
 * @param arr An array of monotonically increasing numbers
 * @param value A number
 * @returns Position in `arr` at which `value` could be inserted while maintaining `arr`'s sort order
 */
export function lowerBound(arr: number[], value: number): number {
    if (arr.length === 0) {
        return 0;
    }
    let low = 0;
    let high = arr.length;

    while (low < high) {
        const mid = Math.floor((low + high) / 2);
        if (arr[mid] < value) {
            low = mid + 1;
        }
        else {
            high = mid;
        }
    }
    return high;
}

/**
 * Piecewise cubic interpolation, using specified derivatives at the interpolating points
 * @param x Independent variable values to interpolate between
 * @param y Dependent variable values to interpolate between
 * @param m Derivatives of the interpolant at the points described by x, y
 * @param xI Independent variable values to interpolate
 * @returns Dependent variable values corresponding to xI, computed using a shape-preserving piecewise cubic Hermite spline
 */
export function piecewiseCubic(x: number[], y: number[], m: number[], xI: number[]): number[] {
    const yI = new Array(xI.length);
    for (let i = 0; i < x.length - 1; ++i) {
        const dX = x[i + 1] - x[i];
        const dY = y[i + 1] - y[i];
        const c = (dY / dX - m[i]) / dX;
        const d = (m[i] + m[i + 1] - (2 * dY) / dX) / (dX * dX);

        const leftIndex = lowerBound(xI, x[i]);
        let rightIndex = lowerBound(xI, x[i + 1]);
        if (i === x.length - 2 && xI[rightIndex] === x[i + 1]) {
            ++rightIndex;
        }
        const xISubset = xI.slice(leftIndex, rightIndex);
        const yISubset = xISubset.map((v) => y[i] + (v - x[i]) * (m[i] + (v - x[i]) * (c + d * (v - x[i + 1]))));
        for (let j = 0; j < yISubset.length; ++j) {
            yI[j + leftIndex] = yISubset[j];
        }
    }
    return yI;
}

/**
 * Similar to Octave's interp1("pchip")
 * @param x Independent variable values to interpolate between
 * @param y Dependent variable values to interpolate between
 * @param xI Independent variable values to interpolate
 * @returns Dependent variable values corresponding to xI, computed using a shape-preserving piecewise cubic Hermite spline
 */
export function pchip(x: number[], y: number[], xI: number[]) {
    return piecewiseCubic(x, y, dpchim(x, y), xI);
}

export class Pchip {
    private deriv: number[];
    constructor(private x: number[], private y: number[]) {
        this.deriv = dpchim(x, y);
    }

    evaluate(xI: number | number[]) {
        if (Array.isArray(xI)) {
            return piecewiseCubic(this.x, this.y, this.deriv, xI);
        }
        else {
            return piecewiseCubic(this.x, this.y, this.deriv, [xI])[0];
        }
    }
}

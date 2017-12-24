
package jsat.linear;

import jsat.math.MathTricks;

import java.util.Arrays;
import java.util.List;

import static java.lang.Math.*;

/**
 * A vector implementation that is dense, meaning all values are allocated -
 * even if their values will be implicitly zero.
 *
 * @author Edward Raff
 */
public class DenseVector extends Vec {

    private static final long serialVersionUID = -889493251793828934L;
    protected final double[] array;


    private int startIndex;
    private int endIndex;


    final static int SMALL_THRESH = 4;

    public static DenseVector a(int length) {
        if (length <= SMALL_THRESH) {
            return new SmallDenseVector(length);
        } else {
            return new DenseVector(length);
        }
    }

    public static DenseVector a(double[] array) {
        if (array.length <= SMALL_THRESH) {
            return new SmallDenseVector(array);
        } else {
            return new DenseVector(array);
        }
    }

    public static class SmallDenseVector extends DenseVector {

        public SmallDenseVector(int size) {
            super(size);
        }

        public SmallDenseVector(double[] array) {
            super(array);
        }

        @Override
        protected void clearStats() {
            //nothing
        }

        @Override
        protected void clearStats(double incoming) {
            //nothing
        }

        @Override
        protected double setStats(double result, int index) {
            //nothing
            return result;
        }
    }

    /**
     * Creates a new Dense Vector of zeros
     *
     * @param length the length of the vector
     */
    private DenseVector(int length) {
        if (length < 0)
            throw new ArithmeticException("You can not have a negative dimension vector");
        array = new double[length];
        startIndex = 0;
        endIndex = array.length;
    }

    /**
     * Creates a new vector of the length of the given list, and values copied
     * over in order.
     *
     * @param list the list of values to copy into a new vector
     */
    public DenseVector(List<Double> list) {
        this.array = new double[list.size()];
        for (int i = 0; i < list.size(); i++)
            this.array[i] = list.get(i);
        startIndex = 0;
        endIndex = this.array.length;
    }

    /**
     * Creates a new Dense Vector that uses the given array as its values. Its
     * values will not be copied, and raw access and mutations tot he given
     * array may occur.
     *
     * @param array the backing array to use for a new vector of the same length
     */
    private DenseVector(double[] array) {
        this(array, 0, array.length);
    }

    /**
     * Creates a new Dense Vector that uses the given array as its values. Its
     * values will not be copied, and raw access and mutations tot he given
     * array may occur.
     *
     * @param array the backing array to use for a new vector
     * @param start the first index in the array, inclusive, to mark the start
     *              of the vector.
     * @param end   the last index in the array, exclusive, to mark the end of the
     *              vector.
     */
    public DenseVector(double[] array, int start, int end) {
        this.array = array;
        this.startIndex = start;
        this.endIndex = end;
    }

    /**
     * Creates a new Dense Vector that contains a copy of the values in the
     * given vector
     *
     * @param toCopy the vector to copy
     */
    public DenseVector(Vec toCopy) {
        this(toCopy.length());
        for (IndexValue iv : toCopy)
            set(iv.index, iv.value);
    }


    @Override
    public int length() {
        return (endIndex - startIndex);
    }

    @Override
    public double get(int index) {
        return array[index + startIndex];
    }

    @Override
    public void set(int index, double val) {
        clearStats(val);
        array[index + startIndex] = val;
    }

    @Override
    public double min() {
        if (stats != null) { double x = stats[MIN]; if (x == x)
            return x; }

        double result = Double.POSITIVE_INFINITY;
        for (int i = startIndex; i < endIndex; i++)
            result = Math.min(result, array[i]);

        return setStats(result, MIN);
    }


    @Override
    public double max() {
        if (stats != null) {
            double x = stats[MAX];
            if (x == x)
                return x;
        }

        double result = Double.NEGATIVE_INFINITY;
        for (int i = startIndex; i < endIndex; i++)
            result = Math.max(result, array[i]);

        return setStats(result, MAX);
    }

    @Override
    public double sum() {
        if (stats != null) {
            double x = stats[SUM];
            if (x == x)
                return x;
        }
        /*
         * Uses Kahan summation algorithm, which is more accurate then
         * naively summing the values in floating point. Though it
         * does not guarenty the best possible accuracy
         *
         * See: http://en.wikipedia.org/wiki/Kahan_summation_algorithm
         */

        //double sum = sumKahan();
        double sum = sumNaive();

        return setStats(sum, SUM);
    }

    private double sumNaive() {
        double sum = 0;
        for (int i = startIndex; i < endIndex; i++) {
            sum += array[i];
        }
        return sum;
    }

    private double sumKahan() {
        double sum = 0;
        double c = 0;
        for (int i = startIndex; i < endIndex; i++) {
            double y = array[i] - c;
            double t = sum + y;
            c = (t - sum) - y;
            sum = t;
        }
        return sum;
    }

    @Override
    public double median() {
        double[] copy = Arrays.copyOfRange(array, startIndex, endIndex);

        Arrays.sort(copy);

        if (copy.length % 2 == 1)
            return copy[copy.length / 2];
        else
            return copy[copy.length / 2] / 2 + copy[copy.length / 2 + 1] / 2;//Divisions by 2 then add is more numericaly stable
    }

    @Override
    public double skewness() {
        double mean = mean();

        double tmp = 0;

        for (int i = startIndex; i < endIndex; i++)
            tmp += MathTricks.cube(array[i] - mean);

        double s1 = tmp / (MathTricks.cube(standardDeviation()) * (array.length - 1));

        if (array.length >= 3)//We can use the bias corrected formula
            return sqrt(array.length * (array.length - 1)) / (array.length - 2) * s1;

        return s1;
    }

    @Override
    public double kurtosis() {
        double mean = mean();

        double tmp = 0;

        for (int i = startIndex; i < endIndex; i++)
            tmp += pow(array[i] - mean, 4);
        tmp /= length();

        return tmp / pow(variance(), 2) - 3;
    }


    @Override
    public DenseVector sortedCopy() {
        double[] copy = Arrays.copyOfRange(array, startIndex, endIndex);

        Arrays.sort(copy);

        return a(copy);
    }

    @Override
    public double variance() {
        if (stats!=null) {
            double x = stats[VAR];
            if (x==x)
                return x;
        }

        double mu = mean();
        double tmp = 0;

        double N = length();


        for (int i = startIndex; i < endIndex; i++)
            tmp += MathTricks.sqr(array[i] - mu) / N;

        return setStats(tmp, VAR);
    }

    @Override
    public double dot(Vec v) {
        if (this.length() != v.length())
            throw new ArithmeticException("Vectors must have the same length");

        if (v.isSparse())
            return v.dot(this);

        double dot = 0;
        for (int i = startIndex; i < endIndex; i++)
            dot += array[i] * v.get(i - startIndex);

        return dot;
    }

    public DenseVector deepCopy() {
        return a(Arrays.copyOf(array, array.length));
    }

    @Override
    public void multiply(double c, Matrix A, Vec b) {

        if (this.length() != A.rows())
            throw new ArithmeticException("Vector x Matrix dimensions do not agree [1," + this.length() + "] x [" + A.rows() + ", " + A.cols() + "]");
        if (b.length() != A.cols())
            throw new ArithmeticException("Destination vector is not the right size");

        for (int i = 0; i < this.length(); i++) {
            double this_i = c * this.array[i + this.startIndex];
            for (int j = 0; j < A.cols(); j++)
                b.increment(j, this_i * A.get(i, j));
        }
    }




    @Override
    public void mutableAdd(double c) {
        if (Math.abs(c) > Double.MIN_NORMAL) {
            clearStats();
            for (int i = startIndex; i < endIndex; i++)
                array[i] += c;
        }
    }

    @Override
    public void mutableAdd(double c, Vec b) {
        if (Math.abs(c) > Double.MIN_NORMAL) {

            if (this.length() != b.length())
                throw new ArithmeticException("Can not add vectors of unequal length");

            clearStats();
            double[] a = this.array;
            if (b instanceof DenseVector) {
                DenseVector db = (DenseVector) b;
                for (int i = startIndex; i < endIndex; i++)
                    a[i] += c * db.array[i];
            } else {
                if (b.isSparse())
                    for (IndexValue iv : b)
                        a[iv.index] += c * iv.value;
                else
                    for (int i = startIndex; i < endIndex; i++)
                        a[i] += c * b.get(i);
            }
        }
    }


    @Override
    public void mutableMultiply(double c) {
        if (Math.abs(c-1) > Double.MIN_NORMAL) {
            clearStats();
            for (int i = startIndex; i < endIndex; i++)
                array[i] *= c;
        }

    }

    @Override
    public void mutableDivide(double c) {
        if (Math.abs(c-1) > Double.MIN_NORMAL) {
            clearStats();
            for (int i = startIndex; i < endIndex; i++)
                array[i] /= c;
        }
    }

    @Override
    public double pNormDist(double p, Vec y) {
        if (this.length() != y.length())
            throw new ArithmeticException("Vectors must be of the same length");

        double norm = 0;
        if (y.isSparse()) {
            int lastIndx = -1;
            for (IndexValue iv : y) {
                for (int i = lastIndx + 1; i < iv.index; i++)//add all the indecies we skipped
                    norm += Math.pow(Math.abs(array[i]), p);
                lastIndx = iv.index;
                //add current
                norm += Math.pow(Math.abs(array[iv.index] - iv.value), p);
            }

            //Tailing zeros
            for (int i = lastIndx + 1; i < y.length(); i++)
                norm += Math.pow(Math.abs(array[i]), p);
        } else {
            for (int i = startIndex; i < endIndex; i++)
                norm += Math.pow(Math.abs(array[i] - y.get(i)), p);
        }
        return Math.pow(norm, 1.0 / p);
    }

    @Override
    public double pNorm(double p) {
        if (p <= 0)
            throw new IllegalArgumentException("norm must be a positive value, not " + p);
        double result = 0;
        if (p == 1) {
            for (int i = startIndex; i < endIndex; i++)
                result += abs(array[i]);
        } else if (p == 2) {
            for (int i = startIndex; i < endIndex; i++)
                result += array[i] * array[i];
            result = Math.sqrt(result);
        } else if (Double.isInfinite(p)) {
            for (int i = startIndex; i < endIndex; i++)
                result = Math.max(result, abs(array[i]));
        } else {
            for (int i = startIndex; i < endIndex; i++)
                result += Math.pow(Math.abs(array[i]), p);
            result = pow(result, 1 / p);
        }
        return result;
    }

    @Override
    public Vec clone() {
        DenseVector copy = a(length());

        System.arraycopy(this.array, startIndex, copy.array, 0, length());

        return copy;
    }

    @Override
    public void normalize() {
        double sum = 0;

        for (int i = startIndex; i < endIndex; i++)
            sum += array[i] * array[i];

        sum = Math.sqrt(sum);

        mutableDivide(Math.max(sum, 1e-10));
    }

    @Override
    public void mutablePairwiseMultiply(Vec b) {
        if (this.length() != b.length())
            throw new ArithmeticException("Vectors must have the same length");
        for (int i = startIndex; i < endIndex; i++)
            this.array[i] *= b.get(i);
    }

    @Override
    public void mutablePairwiseDivide(Vec b) {
        if (this.length() != b.length())
            throw new ArithmeticException("Vectors must have the same length");
        for (int i = startIndex; i < endIndex; i++)
            this.array[i] /= b.get(i);
    }

    @Override
    public boolean equals(Object obj) {
        if (!(obj instanceof Vec))
            return false;
        Vec otherVec = (Vec) obj;

        if (this.length() != otherVec.length())
            return false;
        for (int i = startIndex; i < endIndex; i++)
            if (this.get(i) != otherVec.get(i))
                if (Double.isNaN(this.get(i)) && Double.isNaN(otherVec.get(i)))//NaN != NaN is always true, so check special
                    return true;
                else
                    return false;

        return true;
    }

    @Override
    public boolean equals(Object obj, double range) {
        if (!(obj instanceof Vec))
            return false;
        Vec otherVec = (Vec) obj;
        range = Math.abs(range);

        if (this.length() != otherVec.length())
            return false;
        for (int i = startIndex; i < endIndex; i++)
            if (Math.abs(this.get(i) - otherVec.get(i)) > range)
                if (Double.isNaN(this.get(i)) && Double.isNaN(otherVec.get(i)))//NaN != NaN is always true, so check special
                    return true;
                else
                    return false;

        return true;
    }

    /**
     * Returns a new dense vector backed by the given array. This is a weak
     * reference, the given array should no longer be altered - as it will
     * effect the values of the dense vector.
     *
     * @param array the array to use as the backing of a dense vector
     * @return a Dense Vector that is backed using the given array
     */
    public static DenseVector toDenseVec(double... array) {
        return a(array);
    }

    @Override
    public double[] arrayCopy() {
        return Arrays.copyOfRange(array, startIndex, endIndex);
    }

    @Override
    public boolean isSparse() {
        return false;
    }
}

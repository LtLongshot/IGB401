using System;
using System.Numerics;
using System.Threading.Tasks;


namespace DigitalMusicAnalysis
{
    public class timefreq
    {
        public float[][] timeFreqData;
        public int wSamp;
        public Complex[] twiddles;

        //This is taking forever needs the improvement
        public timefreq(float[] x, int windowSamp)
        {
            //int ii;
            double pi = 3.14159265;
            Complex i = Complex.ImaginaryOne;
            this.wSamp = windowSamp;
            twiddles = new Complex[wSamp];

            Parallel.For(0, wSamp, ii =>
            {
                double a = 2 * pi * ii / (double)wSamp;
                twiddles[ii] = Complex.Pow(Complex.Exp(-i), (float)a);
            });

            timeFreqData = new float[wSamp/2][];

            int nearest = (int)Math.Ceiling((double)x.Length / (double)wSamp);
            nearest = nearest * wSamp;

            Complex[] compX = new Complex[nearest];
            for (int kk = 0; kk < nearest; kk++)
            {
                if (kk < x.Length)
                {
                    compX[kk] = x[kk];
                }
                else
                {
                    compX[kk] = Complex.Zero;
                }
            }


            int cols = 2 * nearest /wSamp;

            for (int jj = 0; jj < wSamp / 2; jj++)
            {
                timeFreqData[jj] = new float[cols];
            }

            timeFreqData = stft(compX, wSamp);
	
        }

        float[][] stft(Complex[] x, int wSamp)
        {
            int N = x.Length;
            float fftMax = 0;

            float[][] Y = new float[wSamp / 2][];

            for (int ll = 0; ll < wSamp / 2; ll++)
            {
                Y[ll] = new float[2 * (int)Math.Floor((double)N / (double)wSamp)];
            }

            //Complex[] temp = new Complex[wSamp];
            Complex[][] tempFFT = new Complex[(int)(2 * Math.Floor((double)N / (double)wSamp) - 1)][];
            Complex[][] tempTemp = new Complex[(int)(2 * Math.Floor((double)N / (double)wSamp) - 1)][];

            for (int ii = 0; ii < 2 * Math.Floor((double)N / (double)wSamp) - 1; ii++)
            {
                Complex[] temp = new Complex[wSamp];
                for (int jj = 0; jj < wSamp; jj++)
                {
                    //outer not safe because of this
                    temp[jj] = x[ii * (wSamp / 2) + jj];
                    }
                tempTemp[ii] = temp;
            }

            Parallel.For(0, (int)(2 * Math.Floor((double)N / (double)wSamp) - 1), ii =>
            {
                tempFFT[ii] = fft(tempTemp[ii]);
            });
            //}

            //gets rid of flow dependency
            //for (int ii = 0; ii < 2 * Math.Floor((double)N / (double)wSamp) - 1; ii++)
            //{
            //paralising this could increase times
            Parallel.For(0, (int)(2 * Math.Floor((double)N / (double)wSamp) - 1), ii => {
                for (int kk = 0; kk < wSamp / 2; kk++) {
                    //Parallel.For(0, wSamp / 2, kk =>
                    //{
                        Y[kk][ii] = (float)Complex.Abs(tempFFT[ii][kk]);

                        if (Y[kk][ii] > fftMax)
                        {
                            fftMax = Y[kk][ii];
                        }
                    //});
                }
            });

            //NOT SAFE could try parallel
            for (int ii = 0; ii < 2 * Math.Floor((double)N / (double)wSamp) - 1; ii++)
            {
                for (int kk = 0; kk < wSamp / 2; kk++)
                {
                    Y[kk][ii] /= fftMax;
                }
            }

            return Y;
        }

        //Complex[] fft(Complex[] x)
        //{
        //    //int ii = 0;
        //    //int kk = 0;
        //    int N = x.Length;

        //    Complex[] Y = new Complex[N];

        //    // NEED TO MEMSET TO ZERO?

        //    if (N == 1)
        //    {
        //        Y[0] = x[0];
        //    }
        //    else
        //    {

        //        Complex[] E = new Complex[N / 2];
        //        Complex[] O = new Complex[N / 2];
        //        Complex[] even = new Complex[N / 2];
        //        Complex[] odd = new Complex[N / 2];

        //        for (int ii = 0; ii < N; ii++)
        //        {
        //            //Parallel.For(0, N, ii =>
        //            //{
        //            if (ii % 2 == 0)
        //            {
        //                even[ii / 2] = x[ii];
        //            }
        //            if (ii % 2 == 1)
        //            {
        //                odd[(ii - 1) / 2] = x[ii];
        //            }
        //            //});

        //        }

        //        E = fft(even);
        //        O = fft(odd);

        //        for (int kk = 0; kk < N; kk++)
        //        {
        //            //Parallel.For(0, N, kk =>
        //            //{
        //            Y[kk] = E[(kk % (N / 2))] + O[(kk % (N / 2))] * twiddles[kk * wSamp / N];

        //            //});
        //        }
        //    }

        //    return Y;
        //}
        public static int BitReverse(int n, int bits)
        {
            int reversedN = n;
            int count = bits - 1;

            n >>= 1;
            while (n > 0)
            {
                reversedN = (reversedN << 1) | (n & 1);
                count--;
                n >>= 1;
            }

            return ((reversedN << count) & ((1 << bits) - 1));
        }

        public Complex[] fft(Complex[] x)
        {
            int length = x.Length;
            Complex[] output = new Complex[length];

            int bits = (int)Math.Log(length, 2);
            for (int j = 0; j < length; j++)
            {
                int swapPos = BitReverse(j, bits);
                output[j] = x[swapPos];
            }
            for (int N = 2; N <= length; N <<= 1)
            {
                for (int i = 0; i < length; i += N)
                {
                    for (int k = 0; k < N / 2; k++)
                    {
                        int evenIndex = i + k;
                        int oddIndex = i + k + (N / 2);
                        var even = output[evenIndex];
                        var odd = output[oddIndex];
                        output[evenIndex] = even + odd * twiddles[k * (wSamp / N)];
                        output[oddIndex] = even + odd * twiddles[(k + (N / 2)) * (wSamp / N)];
                    }
                }
            }

            return output;
        }
    }
    }

using System;
using System.Linq;

namespace vaml_scripts
{
    class Program
    {
        static void Main(string[] args)
        {
            if (args[0].Equals("otsus"))
            {
                string[] histStr = args[1].Split(',');
                OtsusMethod(histStr);
            }
            if (args[0].Equals("gaussian"))
            {
                string[] histStr = args[1].Split(',');
                var modes = 1;;
                if(args[2] != null) modes = Int32.Parse(args[2]);
                Gaussian(histStr,modes);
            }
            if (args[0].Equals("bhdistance")){
                string[] p = args[1].Split(',');
                string[] q = args[2].Split(',');
                double sigma = Double.Parse(args[3]);
                BHDistance(p,q,sigma);
            }
            if(args[0].Equals("normaliseweights")){
                var amount = Int32.Parse(args[1]);
                var weights = new double[amount];
                for(int i=0;i<amount;i++){
                    weights[i] = Double.Parse(args[i+2]);
                }
                NormaliseWeights(weights);
            }
            if(args[0].Equals("gaussianprob")) {
                var toCheck = Double.Parse(args[1]);
                var sigma = Double.Parse(args[2]);
                var u = Double.Parse(args[3]);
                GaussianProb(toCheck,sigma,u);
            }
            if(args[0].Equals("gaussianupdate")){
                var alpha = Double.Parse(args[1]);
                var sigma = Double.Parse(args[2]);
                var u = Double.Parse(args[3]);
                var w = Double.Parse(args[4]);
                var x = Double.Parse(args[5]);
                GaussianUpdate(alpha,sigma,u,w,x);
            }
        }
        static void OtsusMethod(string[] histogramStr)
        {
            var histogram = histogramStr.Select(x => Double.Parse(x)).ToArray();
            var normalisedHist = NormaliseWeights(histogram);
            //run iterations
            double[] sigmas = new double[histogram.Length];
            for (int T = 0; T < histogram.Length; T++)
            {
                //calc weights
                double w1 = 0;
                for (int i = 0; i <= T; i++)
                {
                    w1 += normalisedHist[i];
                }
                double w2 = 1 - w1;
                Console.WriteLine("Iteration: " + T + " w1: " + w1 + " w2: " + w2);
                double initSum = 0;
                for (int i = 0; i <= T; i++)
                {
                    initSum += (i * normalisedHist[i]);
                }
                double u1 = initSum * (1 / w1);
                initSum = 0;
                for (int i = T + 1; i < normalisedHist.Length; i++)
                {
                    initSum += (i * normalisedHist[i]);
                }
                double u2 = initSum * (1 / w2);
                Console.WriteLine("Iteration: " + T + " u1: " + u1 + " u2: " + u2);
                double sigma = w1 * w2 * Math.Pow(u1 - u2, 2);
                Console.WriteLine("Iteration: " + T + " sigma: " + sigma);
                sigmas[T] = sigma;
            }
            int maxIter = Array.IndexOf(sigmas, sigmas.Max());
            Console.WriteLine("T*: " + maxIter + " sigma: " + sigmas.Max());
        }

        static void Gaussian(string[] valuesStr, int modes=1)
        {
            var values = valuesStr.Select(x => Int32.Parse(x)).ToArray();
            var chunkSize = values.Length/modes;
            for(int i=0;i<modes;i++){
                var chunk = values.Skip(chunkSize*i).Take(chunkSize).ToArray();
                double u = (double)chunk.Average();
                foreach(int val in chunk) Console.Write(val+" ");
                Console.WriteLine("mode: "+(i+1)+" u: "+u);
                double sigma = 0;
                for(int j=0;j<chunk.Length;j++){
                    sigma += Math.Pow(chunk[j]-u,2);
                }
                sigma = sigma/(chunk.Length-1);
                sigma = Math.Sqrt(sigma);
                Console.WriteLine("mode: "+(i+1)+" sigma: "+sigma);
            }
        }
        static void GaussianProb(double toCheck, double sigma, double u){
            var prob = 1/Math.Sqrt(2*Math.PI*Math.Pow(sigma,2));
            prob *= Math.Pow(Math.E,-Math.Pow(toCheck-u,2)/(2*sigma*sigma));
            Console.WriteLine("Probability: "+prob);
        }
        static void GaussianUpdate(double alpha, double sigma, double u, double w, double x){
            var newU = (1-alpha)*u+alpha*x;
            Console.WriteLine("newU: "+newU);
            var newSigma = Math.Sqrt((1-alpha)*sigma*sigma+alpha*Math.Pow(x-u,2));
            Console.WriteLine("newSigma: "+newSigma);
            var newW = (1-alpha)*w+alpha;
            Console.WriteLine("newW: "+newW);
        }
        static void BHDistance(string[] pstr, string[] qstr, double sigma){
            var p = pstr.Select(x => Double.Parse(x)).ToArray();
            p = NormaliseWeights(p);
            var q = qstr.Select(x => Double.Parse(x)).ToArray();
            q = NormaliseWeights(q);
            var coeff = 0.0;
            for (int i=0;i<p.Length;i++){
                coeff += Math.Sqrt(p[i]*q[i]);
            }
            Console.WriteLine("ro: " +coeff);
            var dBH = -Math.Log(coeff);
            Console.WriteLine("dBH: "+ dBH);
            var w = 1/(Math.Sqrt(2*Math.PI)*sigma);
            w *= Math.Pow(Math.E,-(dBH*dBH)/(2*sigma*sigma));
            Console.WriteLine("weight: " +w);
        }

        static double[] NormaliseWeights(double[] weights){
            var normalised = new double[weights.Length];
            for(int i=0;i<weights.Length;i++){
                var normalW = weights[i]/weights.Sum();
                Console.WriteLine("w"+i+": "+normalW);
                normalised[i] = normalW;
            }
            return normalised;
        }
    }


}

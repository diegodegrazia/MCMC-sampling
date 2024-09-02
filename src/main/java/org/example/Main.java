package org.example;

import org.example.Plots.MonovariatePlot;
import org.example.Plots.MultivariatePlot;
import org.oristool.math.OmegaBigDecimal;
import org.oristool.math.domain.DBMZone;
import org.oristool.math.expression.*;
import org.oristool.math.function.EXP;
import org.oristool.math.function.GEN;
import org.oristool.math.function.PartitionedGEN;

import javax.sound.midi.Soundbank;
import java.math.BigDecimal;
import java.math.MathContext;
import java.sql.SQLOutput;
import java.util.*;

public class Main {

    private static double bandwidth = 0.1;
    public static void main(String[] args) {

        int numSamples = 10000;

        double c = 100; // constant c for the AR method

        // 0 for the hits, 1 for the misses
        int[] arHitsAndMissesMonovariate = new int[2]; // for monovariate
        int[] arHitsAndMissesMultivariateExp = new int[2]; // for multivariate
        int[] mhHitsAndMissesMonovariate = new int[2]; // for monovariate MH
        int[] mhHitsAndMissesMultivariate = new int[2]; // for multivariate MH

        //variables for the multivariate construction
        Variable x = new Variable("x");
        Variable y = new Variable("y");
        Variable z = new Variable("z");


        //GEN arMonovariateProposal = getARMonovariateProposal(x);
        GEN arMultivariateUniformProposal = getARMultivariateUniformProposal(x,y,z);
        GEN multivariateFunction = getMultivariateFunction(x,y,z);
        //GEN monovariateFunction = getMonovariateFunction(x);
        //PartitionedGEN partitionedGEN = getARMonovariatePiecewiseProposal(x);



        // Sampling with MH and AR
        //List<Double> arPiecewiseSamples = acceptanceRejectionPiecewise(monovariateFunction, numSamples, partitionedGEN, c, arHitsAndMissesMonovariate);
        //double[] mhSamples = metropolisHastings(monovariateFunction, numSamples, mhHitsAndMissesMonovariate);
        //List<Double> arSamples = acceptanceRejection(monovariateFunction, numSamples, arMonovariateProposal, c, arHitsAndMissesMonovariate);
        List<double[]> mhMultiSamles = metropolisHastingsMultivariate(multivariateFunction, numSamples, mhHitsAndMissesMultivariate);
        List<double[]> arMultiSamles = acceptanceRejectionMultivariate(multivariateFunction, numSamples, arMultivariateUniformProposal, c, arHitsAndMissesMonovariate);


        // Plot histograms and scatter plots
        //MonovariatePlot.plotPiecewiseAcceptanceRejectionCharts(arPiecewiseSamples, monovariateFunction, partitionedGEN, c, arHitsAndMissesMonovariate[0], arHitsAndMissesMonovariate[1]);
        //MonovariatePlot.plotCharts(mhSamples, monovariateFunction, mhHitsAndMissesMonovariate[0], mhHitsAndMissesMonovariate[1]);
        //MonovariatePlot.plotAcceptanceRejectionCharts(arSamples, monovariateFunction, arMonovariateProposal, c, arHitsAndMissesMonovariate[0], arHitsAndMissesMonovariate[1]);
        MultivariatePlot.plotMultivariateCharts(mhMultiSamles, mhHitsAndMissesMultivariate[0], mhHitsAndMissesMultivariate[1], "Metropolis-Hastings algorithm");
        MultivariatePlot.plotMultivariateCharts(arMultiSamles, arHitsAndMissesMonovariate[0], arHitsAndMissesMonovariate[1], "Acceptance-Rejection uniform algorithm");

    }


    private static List<Double> acceptanceRejection(GEN targetDistribution, int numSamples, GEN uniformProposal, double c, int [] hitsAndMisses) {
        List<Double> samples = new ArrayList<>();
        Random random = new Random();

        Map<Variable, OmegaBigDecimal> variableMap = new HashMap<>();
        Map<Variable, OmegaBigDecimal> uniformVariableMap = new HashMap<>();
        Variable[] dbmUniformVars = uniformProposal.getDomain().getVariables().toArray(new Variable[uniformProposal.getDomain().getVariables().size()]);
        Variable proposaleVariable = dbmUniformVars[1]; // The uniform proposal is supposed to be a monovariate, uniformVars[0] is the DBM reference variable t*
        Variable[] dbmVars = targetDistribution.getDomain().getVariables().toArray(new Variable[targetDistribution.getDomain().getVariables().size()]);
        Variable targetVariable = dbmVars[1];

        for (int i = 0; i < numSamples; i++) {

            double u = random.nextDouble();
            // uniform proposal sample y
            double y = uniformProposal.getDomainsEFT().doubleValue() + random.nextDouble() * (uniformProposal.getDomainsLFT().doubleValue() - uniformProposal.getDomainsEFT().doubleValue());
            variableMap.put(targetVariable, new OmegaBigDecimal(BigDecimal.valueOf(y)));
            uniformVariableMap.put(proposaleVariable, new OmegaBigDecimal(BigDecimal.valueOf(y)));
            double f_y = targetDistribution.getDensity().evaluate(variableMap).doubleValue();
            double g_y = uniformProposal.getDensity().evaluate(uniformVariableMap).doubleValue(); // uniform density

            if (u <= f_y / (c * g_y)) {
                samples.add(y);
                hitsAndMisses[0]++;
            }
            else
                hitsAndMisses[1]++;
        }
        return samples;
    }

    private static List<Double> acceptanceRejectionWithExponential(GEN targetDistribution, int numSamples, GEN expProposal, double c, int[] hitsAndMisses) {
        List<Double> samples = new ArrayList<>();
        Random random = new Random();
        Map<Variable, OmegaBigDecimal> variableMap = new HashMap<>();
        Variable[] dbmVariable = targetDistribution.getDomain().getVariables().toArray(new Variable[targetDistribution.getDomain().getVariables().size()]); // Assume single variable
        Variable variable = dbmVariable[1];

        for (int i = 0; i < numSamples; i++) {
            double u = random.nextDouble();
            // Exponential proposal sample y
            double y = -Math.log(1 - random.nextDouble()) / expProposal.getDensity().getExponentialRate().doubleValue();

            variableMap.put(variable, new OmegaBigDecimal(BigDecimal.valueOf(y)));
            double f_y = targetDistribution.getDensity().evaluate(variableMap).doubleValue();
            double g_y = expProposal.getDensity().evaluate(variableMap).doubleValue(); // Exponential density

            if (u <= f_y / (c * g_y)) {
                samples.add(y);
                hitsAndMisses[0]++;
            } else {
                hitsAndMisses[1]++;
            }
        }
        return samples;
    }

    private static List<double[]> acceptanceRejectionMultivariate(GEN targetDistribution, int numSamples, GEN uniformProposal, double c, int[] hitsAndMisses) {
        List<double[]> samples = new ArrayList<>();
        Random random = new Random();

        Map<Variable, OmegaBigDecimal> variableMap = new HashMap<>();

        // Get variables from the target and uniform proposal distributions
        Variable[] vars = targetDistribution.getDomain().getVariables().toArray(new Variable[targetDistribution.getDomain().getVariables().size() - 1]);
        List<Variable> variables = new ArrayList<>(Arrays.asList(vars));
        Variable dbmReferenceVariable = variables.remove(0);

        for (int i = 0; i < numSamples; i++) {
            double u = random.nextDouble();
            double[] y = new double[variables.size()];  // Proposal sample in each dimension

            // Generate a uniform proposal sample for each dimension
            for (int d = 0; d < variables.size(); d++) {
                y[d] = uniformProposal.getDomain().getCoefficient(dbmReferenceVariable, variables.get(d)).negate().doubleValue() + random.nextDouble() *
                        (uniformProposal.getDomain().getCoefficient(variables.get(d), dbmReferenceVariable).doubleValue() - uniformProposal.getDomain().getCoefficient(dbmReferenceVariable, variables.get(d)).negate().doubleValue());
                variableMap.put(variables.get(d), new OmegaBigDecimal(BigDecimal.valueOf(y[d])));
            }

            //check if the sample is in the DBM domain of the target, if not discard it immediately
            if(targetDistribution.getDomain().contains(variableMap)) {

                // Evaluate the densities
                double targetEvaluation = targetDistribution.getDensity().evaluate(variableMap).doubleValue();
                double proposalEvaluation = uniformProposal.getDensity().evaluate(variableMap).doubleValue(); // Uniform density

                // Acceptance-Rejection criterion
                if (u <= targetEvaluation / (c * proposalEvaluation)) {
                    samples.add(y);
                    hitsAndMisses[0]++;  // Increment hits
                } else {
                    hitsAndMisses[1]++;  // Increment misses
                }
            }
            else {
                hitsAndMisses[1]++;  // Increment misses
            }
        }
        return samples;
    }

    private static List<double[]> acceptanceRejectionMultivariateWithExponential(GEN targetDistribution, int numSamples, GEN expProposal, double c, int[] hitsAndMissesExp) {
        List<double[]> samples = new ArrayList<>();
        Random random = new Random();

        Map<Variable, OmegaBigDecimal> variableMap = new HashMap<>();
        Variable[] vars = targetDistribution.getDomain().getVariables().toArray(new Variable[targetDistribution.getDomain().getVariables().size() - 1]);
        List<Variable> variables = new ArrayList<>(Arrays.asList(vars));
        variables.remove(0);

        for (int i = 0; i < numSamples; i++) {
            double u = random.nextDouble();
            double[] y = new double[variables.size()];  // Proposal sample in each dimension

            // Generate an exponential proposal sample for each dimension
            for (int d = 0; d < variables.size(); d++) {
                ExponentialTerm et = (ExponentialTerm) expProposal.getDensity().getExmonomials().get(d).getAtomicTerms().get(0);
                y[d] = -Math.log(1 - random.nextDouble()) / et.getLambda().doubleValue();
                variableMap.put(variables.get(d), new OmegaBigDecimal(BigDecimal.valueOf(y[d])));
            }

            // Evaluate the densities
            double targetEvaluation = targetDistribution.getDensity().evaluate(variableMap).doubleValue();
            double proposalEvaluation = expProposal.getDensity().evaluate(variableMap).doubleValue();

            // Acceptance-Rejection criterion
            if (u <= targetEvaluation / (c * proposalEvaluation)) {
                samples.add(y);
                hitsAndMissesExp[0]++;  // Increment hits
            } else {
                hitsAndMissesExp[1]++;  // Increment misses
            }
        }
        return samples;
    }

    private static double[] metropolisHastings(GEN targetDistribution, int numSamples, int[] mhHitsAndMisses) {
        double[] samples = new double[numSamples];
        Random random = new Random();
        MathContext mc = MathContext.DECIMAL128;
        OmegaBigDecimal current = new OmegaBigDecimal(new BigDecimal(targetDistribution.getDomainsEFT().doubleValue() + random.nextDouble(), mc));
        OmegaBigDecimal proposal;
        OmegaBigDecimal alpha;
        Variable[] vars = targetDistribution.getDomain().getVariables().toArray(new Variable[targetDistribution.getDomain().getVariables().size()]);
        Variable targetVariable = vars[1];
        double stdDeviation = 1; // Bigger values for wider boundaries of the variable, otherwise I get "mountains like" histograms.

        Map<Variable, OmegaBigDecimal> variableMap = new HashMap<>();

        for (int i = 0; i < numSamples; i++) {

            proposal = new OmegaBigDecimal(BigDecimal.valueOf(random.nextGaussian()*stdDeviation + current.doubleValue()));

            if (proposal.compareTo(targetDistribution.getDomainsEFT()) < 0  || proposal.compareTo(targetDistribution.getDomainsLFT()) > 0) { // if proposal < minBound or > maxBound then I discard it and stay in the current state.
                samples[i] = current.doubleValue();
                mhHitsAndMisses[1]++;
                continue;
            }

            variableMap.put(targetVariable, current);
            OmegaBigDecimal currentDensityValue = targetDistribution.getDensity().evaluate(variableMap);

            variableMap.put(targetVariable, proposal);
            OmegaBigDecimal proposalDensityValue = targetDistribution.getDensity().evaluate(variableMap);

            alpha = proposalDensityValue.divide(currentDensityValue.bigDecimalValue(), mc).min(OmegaBigDecimal.ONE);

            if (new BigDecimal(random.nextDouble(), mc).compareTo(alpha.bigDecimalValue()) < 0) { // Accept the proposal state
                mhHitsAndMisses[0]++;
                current = proposal;
            }
            else { // Refuse the proposal state
                mhHitsAndMisses[1]++;
            }
            samples[i] = current.doubleValue();
        }
        return samples;
    }

    private static List<double[]> metropolisHastingsMultivariate(GEN multivariate, int numSamples, int[] mhMultiHitsAndMisses) {
        Variable[] vars = multivariate.getDomain().getVariables().toArray(new Variable[multivariate.getDomain().getVariables().size() - 1]);
        List<Variable> variables = new ArrayList<>(Arrays.asList(vars));
        Variable dbmReferenceVariable = variables.remove(0); // t* is its name if printed
        Random random = new Random();
        List<double[]> samples = new ArrayList<>();
        MathContext mc = MathContext.DECIMAL128;
        OmegaBigDecimal[] currents = new OmegaBigDecimal[variables.size()];

        double[] stdDeviation = new double[variables.size()];
        stdDeviation[0] = 0.2;
        stdDeviation[1] = 0.2;
        stdDeviation[2] = 0.2;

        for (int i = 0; i < variables.size(); i++) {
            currents[i] = new OmegaBigDecimal(new BigDecimal(multivariate.getDomain().getCoefficient(dbmReferenceVariable, variables.get(i)).negate().doubleValue() + random.nextDouble() * (multivariate.getDomain().getCoefficient(variables.get(i), dbmReferenceVariable).doubleValue() - multivariate.getDomain().getCoefficient(dbmReferenceVariable, variables.get(i)).negate().doubleValue()), mc));
        }

        OmegaBigDecimal[] proposals = new OmegaBigDecimal[variables.size()];
        OmegaBigDecimal alpha;

        Map<Variable, OmegaBigDecimal> currentVariableMap = new HashMap<>(); // for current density evaluation
        Map<Variable, OmegaBigDecimal> proposalVariableMap = new HashMap<>(); // for proposal density evaluation

        mhMultiHitsAndMisses[0] = 0; // hits
        mhMultiHitsAndMisses[1] = 0; // misses

        outerLoop:
        for (int i = 0; i < numSamples; i++) {
            // Generate proposal by adding a small random perturbation
            for (int d = 0; d < variables.size(); d++) {
                proposals[d] = new OmegaBigDecimal(BigDecimal.valueOf(random.nextGaussian() * stdDeviation[d] + currents[d].doubleValue()));
                // Ensure proposal is between the variable bounds
                if (proposals[d].compareTo(multivariate.getDomain().getCoefficient(dbmReferenceVariable, variables.get(d)).negate()) < 0 || proposals[d].compareTo(multivariate.getDomain().getCoefficient(variables.get(d), dbmReferenceVariable)) > 0) {
                    mhMultiHitsAndMisses[1]++;
                    continue outerLoop;
                }
                // Preparing the variable map for current and proposal density evaluation
                currentVariableMap.put(variables.get(d), currents[d]);
                proposalVariableMap.put(variables.get(d), proposals[d]);
            }

            if (multivariate.getDomain().contains(proposalVariableMap)) {
                OmegaBigDecimal currentDensityValue = multivariate.getDensity().evaluate(currentVariableMap);
                OmegaBigDecimal proposalDensityValue = multivariate.getDensity().evaluate(proposalVariableMap);

                // Calculate acceptance ratio alpha
                alpha = proposalDensityValue.divide(currentDensityValue.bigDecimalValue(), mc).min(OmegaBigDecimal.ONE);

                // Accept or reject proposal
                if (new BigDecimal(random.nextDouble(), mc).compareTo(alpha.bigDecimalValue()) < 0) {
                    double[] acceptedSample = new double[variables.size()];
                    for (int j = 0; j < variables.size(); j++) {
                        currents[j] = proposals[j];
                        acceptedSample[j] = currents[j].doubleValue();
                    }
                    samples.add(acceptedSample);
                    mhMultiHitsAndMisses[0]++;
                } else {
                    mhMultiHitsAndMisses[1]++;
                }
            }
            else {
                mhMultiHitsAndMisses[1]++;
            }
        }
        return samples;
    }

    private static double[] generateExpCDFInversionSamples(EXP exp, int numSamples) {
        double[] samples = new double[numSamples];
        Random random = new Random();
        MathContext mc = MathContext.DECIMAL128;
        double uMin = 1 - Math.exp(-exp.getLambda().doubleValue() * exp.getDomainsEFT().doubleValue());
        double uMax = 1 - Math.exp(-exp.getLambda().doubleValue() * exp.getDomainsLFT().doubleValue());
        double u;
        for (int i = 0; i < numSamples; i++) {
            u = uMin + random.nextDouble() * (uMax - uMin);
            samples[i] = -Math.log(1 - u) / exp.getLambda().doubleValue();
        }
        return samples;
    }

    private static double[] generateSumOfExponentialsSamples(EXP exp, int k, int numSamples) {
        double[] samples = new double[numSamples];
        Random random = new Random();

        for (int i = 0; i < numSamples; i++) {
            double sum = 0.0;
            for (int j = 0; j < k; j++) {
                double u = random.nextDouble();
                double expSample = -Math.log(1 - u) / exp.getLambda().doubleValue();
                sum += expSample;
            }
            samples[i] = sum;
        }
        return samples;
    }

    private static GEN getMonovariateFunction(Variable var) {
        //k=4, lambda=12
        Expolynomial expolynomial = new Expolynomial();
        Exmonomial exmonomial = new Exmonomial(BigDecimal.valueOf(3456));
        exmonomial.addAtomicTerm(new MonomialTerm(var, 3));
        Exmonomial exmonomial2 = new Exmonomial(BigDecimal.valueOf(1));
        exmonomial2.addAtomicTerm(new ExponentialTerm(var, BigDecimal.valueOf(12)));
        exmonomial.multiply(exmonomial2);
        expolynomial.addExmonomial(exmonomial);
        GEN erlang = new GEN(new DBMZone(var), expolynomial);
        erlang.getDomain().imposeBound(Variable.TSTAR, var, OmegaBigDecimal.ZERO);
        erlang.getDomain().imposeBound(var, Variable.TSTAR, OmegaBigDecimal.POSITIVE_INFINITY);
        return erlang;
    }

    private static GEN getMultivariateFunction(Variable x, Variable y, Variable z) {

        //Multivariate target as product
        Exmonomial ex = new Exmonomial(BigDecimal.valueOf(1));
        ex.addAtomicTerm(new MonomialTerm(x,1));
        Exmonomial ex2 = new Exmonomial(BigDecimal.valueOf(1));
        ExponentialTerm exp1 = new ExponentialTerm(x, BigDecimal.valueOf(3));
        Exmonomial ex3 = new Exmonomial(BigDecimal.valueOf(1));
        ExponentialTerm exp2 = new ExponentialTerm(y, BigDecimal.valueOf(-2));
        ex2.addAtomicTerm(exp1);
        ex3.addAtomicTerm(exp2);
        ex.multiply(ex2);
        ex.multiply(ex3); // x*e^(-3x+2y)

        Exmonomial ex4 = new Exmonomial(BigDecimal.valueOf(1));
        MonomialTerm mt1 = new MonomialTerm(x, 2);
        ex4.addAtomicTerm(mt1);
        Exmonomial ex5 = new Exmonomial(BigDecimal.valueOf(1));
        MonomialTerm mt2 = new MonomialTerm(y, 2);
        ex5.addAtomicTerm(mt2);
        ex4.multiply(ex5);
        Exmonomial ex6 = new Exmonomial(BigDecimal.valueOf(1));
        ExponentialTerm exp3 = new ExponentialTerm(z, BigDecimal.valueOf(5));
        ex6.addAtomicTerm(exp3);
        ex4.multiply(ex6);
        Exmonomial ex7 = new Exmonomial(BigDecimal.valueOf(1));
        ExponentialTerm exp4 = new ExponentialTerm(x, BigDecimal.valueOf(3));
        ex7.addAtomicTerm(exp4);
        ex4.multiply(ex7);
        Exmonomial ex8 = new Exmonomial(BigDecimal.valueOf(1));
        ExponentialTerm exp5 = new ExponentialTerm(y, BigDecimal.valueOf(-2));
        ex8.addAtomicTerm(exp5);
        ex4.multiply(ex8);

        Expolynomial expolynomial = new Expolynomial();
        expolynomial.addExmonomial(ex);
        expolynomial.addExmonomial(ex4);

        GEN multivariateFunction = new GEN(new DBMZone(x,y,z), expolynomial);
        multivariateFunction.getDomain().imposeBound(Variable.TSTAR,x, OmegaBigDecimal.ZERO);
        multivariateFunction.getDomain().imposeBound(Variable.TSTAR,y, OmegaBigDecimal.ZERO);
        multivariateFunction.getDomain().imposeBound(Variable.TSTAR,z, OmegaBigDecimal.ZERO);
        multivariateFunction.getDomain().imposeBound(x,Variable.TSTAR, new OmegaBigDecimal(BigDecimal.valueOf(10)));
        multivariateFunction.getDomain().imposeBound(y,Variable.TSTAR, new OmegaBigDecimal(BigDecimal.valueOf(10)));
        multivariateFunction.getDomain().imposeBound(z,Variable.TSTAR, new OmegaBigDecimal(BigDecimal.valueOf(8)));
        System.out.println(multivariateFunction);
        return multivariateFunction;
    }

    private static GEN getARMonovariateProposal(Variable var) {

        return GEN.newUniform(OmegaBigDecimal.ZERO, new OmegaBigDecimal(BigDecimal.valueOf(1.5)));
    }

    private static GEN getARMultivariateUniformProposal(Variable x, Variable y, Variable z) {
        //Multivariate uniform proposal construction for the AR method
        Expolynomial expolynomialUniformProposalFunction = new Expolynomial();

        Exmonomial exUniformProposal = new Exmonomial(BigDecimal.valueOf(0.1));
        exUniformProposal.addAtomicTerm(new MonomialTerm(x, 0));
        Exmonomial ex2UniformProposal = new Exmonomial(BigDecimal.valueOf(0.1));
        ex2UniformProposal.addAtomicTerm(new MonomialTerm(y, 0));
        Exmonomial ex3UniformProposal = new Exmonomial(BigDecimal.valueOf(0.125));
        ex3UniformProposal.addAtomicTerm(new MonomialTerm(z, 0));
        exUniformProposal.multiply(ex2UniformProposal);
        exUniformProposal.multiply(ex3UniformProposal);
        expolynomialUniformProposalFunction.addExmonomial(exUniformProposal);

        GEN arMultivariateUniformProposal = new GEN(new DBMZone(x,y,z), expolynomialUniformProposalFunction);
        arMultivariateUniformProposal.getDomain().imposeBound(x, Variable.TSTAR, new OmegaBigDecimal(BigDecimal.valueOf(10))); // upper bound
        arMultivariateUniformProposal.getDomain().imposeBound(Variable.TSTAR, x, new OmegaBigDecimal(BigDecimal.valueOf(0)).negate()); // lower bound
        arMultivariateUniformProposal.getDomain().imposeBound(y, Variable.TSTAR, new OmegaBigDecimal(BigDecimal.valueOf(10))); // upper bound
        arMultivariateUniformProposal.getDomain().imposeBound(Variable.TSTAR, y, new OmegaBigDecimal(BigDecimal.valueOf(0)).negate()); // lower bound
        arMultivariateUniformProposal.getDomain().imposeBound(z, Variable.TSTAR, new OmegaBigDecimal(BigDecimal.valueOf(8))); // upper bound
        arMultivariateUniformProposal.getDomain().imposeBound(Variable.TSTAR, z, new OmegaBigDecimal(BigDecimal.valueOf(0)).negate()); // lower bound
        return arMultivariateUniformProposal;
    }

    private static double klDivergence(GEN p, GEN q, List<Double> samples) { // Note klDivergence is not symmetric!
        double divergence = 0;
        Variable[] dbmVars = p.getDomain().getVariables().toArray(new Variable[p.getDomain().getVariables().size()]);
        Variable var = dbmVars[1]; // the first is t*
        for (int i=0; i<samples.size(); i++) {
            Map<Variable, OmegaBigDecimal> variableMap = new HashMap<>();
            variableMap.put(var, new OmegaBigDecimal(BigDecimal.valueOf(samples.get(i))));
            divergence += p.getDensity().evaluate(variableMap).doubleValue() * (Math.log(p.getDensity().evaluate(variableMap).doubleValue()/q.getDensity().evaluate(variableMap).doubleValue())/Math.log(2));
        }
        return divergence;
    }

    private static double jsDivergence(GEN p, GEN q, List<Double> samples) {
        double divergence = 0;
        Variable[] dbmVars = p.getDomain().getVariables().toArray(new Variable[p.getDomain().getVariables().size()]);
        Variable var = dbmVars[1]; // the first is t*

        for (int i = 0; i < samples.size(); i++) {
            Map<Variable, OmegaBigDecimal> variableMap = new HashMap<>();
            variableMap.put(var, new OmegaBigDecimal(BigDecimal.valueOf(samples.get(i))));

            double pDensity = p.getDensity().evaluate(variableMap).doubleValue();
            double qDensity = q.getDensity().evaluate(variableMap).doubleValue();
            double mDensity = (pDensity + qDensity) / 2.0;

            // KL divergence part for p with respect to m
            double klPm = pDensity * (Math.log(pDensity / mDensity) / Math.log(2));

            // KL divergence part for q with respect to m
            double klQm = qDensity * (Math.log(qDensity / mDensity) / Math.log(2));

            // Add both parts to the total divergence
            divergence += 0.5 * klPm + 0.5 * klQm;
        }

        return divergence;
    }


    private static List<Double> acceptanceRejectionPiecewise(GEN targetDistribution, int numSamples, PartitionedGEN piecewiseProposal, double c, int[] hitsAndMisses) {
        List<Double> samples = new ArrayList<>();
        Random random = new Random();
        Map<Variable, OmegaBigDecimal> variableMap = new HashMap<>();

        // Split the proposal into its components (Uniform and Exponential)
        GEN uniformProposal = piecewiseProposal.getFunctions().get(0);
        GEN expProposal = piecewiseProposal.getFunctions().get(1);

        // Variables in the target distribution
        Variable[] dbmTargetVariable = targetDistribution.getDomain().getVariables().toArray(new Variable[uniformProposal.getDomain().getVariables().size()]);
        Variable targetVariable = dbmTargetVariable[1]; // the first one is t*

        for (int i = 0; i < numSamples; i++) {
            double u = random.nextDouble();
            double y;
            double g_y;
            double ratio = uniformProposal.integrateOverDomain().doubleValue()/(uniformProposal.integrateOverDomain().doubleValue()+expProposal.integrateOverDomain().doubleValue());

            // Decide which part of the piecewise proposal to sample from
            if (u < ratio) {
                // Sample from Uniform distribution
                y = uniformProposal.getDomainsEFT().doubleValue() + random.nextDouble() * (uniformProposal.getDomainsLFT().doubleValue() - uniformProposal.getDomainsEFT().doubleValue());
                variableMap.put(targetVariable, new OmegaBigDecimal(BigDecimal.valueOf(y)));
                g_y = uniformProposal.getDensity().evaluate(variableMap).doubleValue();
            } else {
                // Sample from Exponential distribution
                double uMin = 1 - Math.exp(-expProposal.getDensity().getExponentialRate().doubleValue() * expProposal.getDomainsEFT().doubleValue());
                double uMax = 1 - Math.exp(-expProposal.getDensity().getExponentialRate().doubleValue() * expProposal.getDomainsLFT().doubleValue());
                y = -Math.log(1 - (uMin + random.nextDouble() * (uMax - uMin))) / expProposal.getDensity().getExponentialRate().doubleValue();
                variableMap.put(targetVariable, new OmegaBigDecimal(BigDecimal.valueOf(y)));
                g_y = expProposal.getDensity().evaluate(variableMap).doubleValue();
            }

            // Evaluate the target distribution
            double f_y = targetDistribution.getDensity().evaluate(variableMap).doubleValue();

            // Apply the Acceptance-Rejection criterion
            if (random.nextDouble() <= f_y / (c * g_y)) {
                samples.add(y);
                hitsAndMisses[0]++; // count as a hit
            } else {
                hitsAndMisses[1]++; // count as a miss
            }
        }
        return samples;
    }

    private static PartitionedGEN getARMonovariatePiecewiseProposal(Variable var) {

        OmegaBigDecimal epsilon = new OmegaBigDecimal(BigDecimal.valueOf(0.05));
        GEN monovariatePiecewise = GEN.newUniform(OmegaBigDecimal.ZERO, new OmegaBigDecimal(BigDecimal.valueOf(0.25)).add(epsilon));
        GEN monovariatePiecewise2 = GEN.newTruncatedExp(var, BigDecimal.valueOf(1/monovariatePiecewise.getDomainsLFT().doubleValue()), monovariatePiecewise.getDomainsLFT(), OmegaBigDecimal.POSITIVE_INFINITY);
        monovariatePiecewise2.getDensity().multiply(BigDecimal.valueOf(Math.exp(1))); //shift it to the right by epsilon.
        List<GEN> monovariatePiecewiseElements = new ArrayList<>();
        monovariatePiecewiseElements.add(monovariatePiecewise);
        monovariatePiecewiseElements.add(monovariatePiecewise2);
        PartitionedGEN monovariateProposal = new PartitionedGEN(monovariatePiecewiseElements);
        return monovariateProposal;
    }

    /*public static GEN kdeEstimation(List<Double> samples, Variable var) {
        double c = 1/(Math.sqrt(2*Math.PI)*samples.size()*bandwidth);
        ExponentialTerm exponentialTerm1 = new ExponentialTerm(var., BigDecimal.valueOf(Math.sqrt(2*bandwidth*bandwidth)));

        for (int i=0; i<samples.size(); i++) {

        }

    }

    */
}



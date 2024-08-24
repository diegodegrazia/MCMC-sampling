package org.example;

import com.google.errorprone.annotations.Var;
import org.example.Plots.MonovariatePlot;
import org.example.Plots.MultivariatePlot;
import org.oristool.math.OmegaBigDecimal;
import org.oristool.math.domain.DBMZone;
import org.oristool.math.expression.*;
import org.oristool.math.function.EXP;
import org.oristool.math.function.GEN;

import javax.swing.*;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.util.*;

public class Main {
    public static void main(String[] args) {

        int numSamples = 10000;

        double c = 24; // constant c for the AR method

        // 0 for the hits, 1 for the misses
        int[] arHitsAndMissesMonovariateUniform = new int[2]; // for monovariate with uniform proposal
        int[] arHitsAndMissesMultivariateExp = new int[2]; // for multivariate with EXP proposal
        int[] mhHitsAndMissesMonovariate = new int[2]; // for monovariate MH
        int[] mhHitsAndMissesMultivariate = new int[2]; // for multivariate MH

        //variables for the multivariate construction
        Variable x = new Variable("x");
        Variable y = new Variable("y");
        Variable z = new Variable("z");


        GEN arMonovariateUniformProposal = getARMonovariateUniformProposal(x);
        //GEN arMultivariateUniformProposal = getARMultivariateUniformProposal(x,y,z);
        //GEN arMultivariateExpProposal = getARMultivariateExpProposal(x,y,z);
        //GEN multivariateFunction = getMultivariateFunction(x,y,z);
        GEN monovariateFunction = getMonovariateFunction(x);


        // Sampling with MH and AR
        double[] mhSamples = metropolisHastings(monovariateFunction, numSamples, mhHitsAndMissesMonovariate);
        //List<Double> arSamples = acceptanceRejection(monovariateFunction, numSamples, arMonovariateUniformProposal, c, arHitsAndMissesMonovariateUniform);
        //List<double[]> mhMultiSamles = metropolisHastingsMultivariate(multivariateFunction, numSamples, mhHitsAndMissesMultivariate);
        //List<double[]> arMultiSamles = acceptanceRejectionMultivariate(multivariateFunction, numSamples, arMultivariateUniformProposal, c, arHitsAndMissesMonovariateUniform);
        //List<double[]> arMultiSamples2 = acceptanceRejectionMultivariateWithExponential(multivariateFunction, numSamples, arMultivariateExpProposal, c, arHitsAndMissesMultivariateExp);


        // Plot histograms and scatter plots.
        MonovariatePlot.plotCharts(mhSamples, monovariateFunction, mhHitsAndMissesMonovariate[0], mhHitsAndMissesMonovariate[1]);
        //MonovariatePlot.plotAcceptanceRejectionCharts(arSamples, monovariateFunction, arMonovariateUniformProposal, c, arHitsAndMissesMonovariateUniform[0], arHitsAndMissesMonovariateUniform[1]);
        //MultivariatePlot.plotMultivariateCharts(mhMultiSamles, mhHitsAndMissesMultivariate[0], mhHitsAndMissesMultivariate[1], "Metropolis-Hastings algorithm");
        //MultivariatePlot.plotMultivariateCharts(arMultiSamles, arHitsAndMissesMonovariateUniform[0], arHitsAndMissesMonovariateUniform[1], "Acceptance-Rejection uniform algorithm");
        //MultivariatePlot.plotMultivariateCharts(arMultiSamples2, arHitsAndMissesMultivariateExp[0], arHitsAndMissesMultivariateExp[1], "Acceptance-Rejection exp algorithm");

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

    private static List<Double> acceptanceRejectionWithExponential(GEN targetDistribution, int numSamples, EXP expProposal, double c, int[] hitsAndMisses) {
        List<Double> samples = new ArrayList<>();
        Random random = new Random();
        Map<Variable, OmegaBigDecimal> variableMap = new HashMap<>();
        Variable targetVariable = targetDistribution.getDomain().getVariables().iterator().next(); // Assume single variable

        for (int i = 0; i < numSamples; i++) {
            double u = random.nextDouble();
            // Exponential proposal sample y
            double y = -Math.log(1 - random.nextDouble()) / expProposal.getLambda().doubleValue();

            variableMap.put(targetVariable, new OmegaBigDecimal(BigDecimal.valueOf(y)));
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
        double stdDeviation = 0.006; // Bigger values for wider boundaries of the variable, otherwise I get "mountains like" histograms.

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

        // TODO: different std deviations for each dimension (makes sense if you think about it).
        double[] stdDeviation = new double[variables.size()];
        stdDeviation[0] = 0.2;
        stdDeviation[1] = 0.2;
        stdDeviation[2] = 0.2;

        for (int i = 0; i < variables.size(); i++) {
            currents[i] = new OmegaBigDecimal(new BigDecimal(multivariate.getDomain().getCoefficient(dbmReferenceVariable, variables.get(i)).negate().doubleValue() + random.nextDouble() * (multivariate.getDomain().getCoefficient(variables.get(i), dbmReferenceVariable).doubleValue() - multivariate.getDomain().getCoefficient(dbmReferenceVariable, variables.get(i)).negate().doubleValue()), mc));
        }

        OmegaBigDecimal[] proposals = new OmegaBigDecimal[variables.size()];
        OmegaBigDecimal alpha;

        Map<Variable, OmegaBigDecimal>[] variableMaps = new HashMap[2]; // for current and proposal
        variableMaps[0] = new HashMap<>(); // 0 for current
        variableMaps[1] = new HashMap<>(); // 1 for proposal

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
                variableMaps[0].put(variables.get(d), currents[d]);
                variableMaps[1].put(variables.get(d), proposals[d]);
            }
            OmegaBigDecimal currentDensityValue = multivariate.getDensity().evaluate(variableMaps[0]);
            OmegaBigDecimal proposalDensityValue = multivariate.getDensity().evaluate(variableMaps[1]);

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

        GEN exp = GEN.newTruncatedExp(var, BigDecimal.valueOf(3), OmegaBigDecimal.ZERO, new OmegaBigDecimal(BigDecimal.valueOf(50)));
        return exp;
    }

    private static GEN getMultivariateFunction(Variable x, Variable y, Variable z) {

        //Multivariate target as product
        Exmonomial ex = new Exmonomial(BigDecimal.valueOf(3.0));
        Exmonomial ex2 = new Exmonomial(BigDecimal.valueOf(1.0));
        Exmonomial ex3 = new Exmonomial(BigDecimal.valueOf(1.0));
        ex.addAtomicTerm(new ExponentialTerm(x, BigDecimal.valueOf(1)));
        ex2.addAtomicTerm(new ExponentialTerm(y, BigDecimal.valueOf(1)));
        ex3.addAtomicTerm(new MonomialTerm(z, 2));
        ex.multiply(ex2);
        ex.multiply(ex3);
        Expolynomial expolynomialFunction = new Expolynomial();
        expolynomialFunction.addExmonomial(ex);
        GEN multivariateFunction = new GEN(new DBMZone(x,y,z), expolynomialFunction);
        multivariateFunction.getDomain().imposeBound(Variable.TSTAR, x, new OmegaBigDecimal(BigDecimal.valueOf(0).negate()));
        multivariateFunction.getDomain().imposeBound(Variable.TSTAR, y, new OmegaBigDecimal(BigDecimal.valueOf(0).negate()));
        multivariateFunction.getDomain().imposeBound(Variable.TSTAR, z, new OmegaBigDecimal(BigDecimal.valueOf(0).negate()));
        multivariateFunction.getDomain().imposeBound(x, Variable.TSTAR, new OmegaBigDecimal(BigDecimal.valueOf(5)));
        multivariateFunction.getDomain().imposeBound(y, Variable.TSTAR, new OmegaBigDecimal(BigDecimal.valueOf(5)));
        multivariateFunction.getDomain().imposeBound(z, Variable.TSTAR, new OmegaBigDecimal(BigDecimal.valueOf(10)));
        System.out.println(multivariateFunction);
        return multivariateFunction;
    }

    private static GEN getARMonovariateUniformProposal(Variable var) {
        GEN monovariateUniformProposal = GEN.newUniform(new OmegaBigDecimal(BigDecimal.valueOf(0)), new OmegaBigDecimal(BigDecimal.valueOf(8)));
        // by default newUniform() sets the variable as x, in this way I impose the requested one.
        Variable[] dbmVars = monovariateUniformProposal.getDomain().getVariables().toArray(new Variable[monovariateUniformProposal.getDomain().getVariables().size()]);
        Variable oldVariable = dbmVars[1];
        monovariateUniformProposal.substitute(oldVariable, var);
        return monovariateUniformProposal;
    }

    private static GEN getARMultivariateExpProposal(Variable x, Variable y, Variable z) {
        //Multivariate exp proposal construction for the AR method
        Expolynomial expolynomialExpProposalFunction = new Expolynomial();

        Exmonomial exExpProposal = new Exmonomial(BigDecimal.valueOf(3));
        exExpProposal.addAtomicTerm(new ExponentialTerm(z, BigDecimal.valueOf(2)));
        expolynomialExpProposalFunction.addExmonomial(exExpProposal);
        GEN arMultivariateExpProposal = new GEN(new DBMZone(x,y,z), expolynomialExpProposalFunction);
        arMultivariateExpProposal.getDomain().imposeBound(Variable.TSTAR, x, new OmegaBigDecimal(BigDecimal.valueOf(0)));
        arMultivariateExpProposal.getDomain().imposeBound(Variable.TSTAR, y, new OmegaBigDecimal(BigDecimal.valueOf(0)));
        arMultivariateExpProposal.getDomain().imposeBound(Variable.TSTAR, z, new OmegaBigDecimal(BigDecimal.valueOf(0)));
        arMultivariateExpProposal.getDomain().imposeBound(x, Variable.TSTAR, new OmegaBigDecimal(BigDecimal.valueOf(50)));
        arMultivariateExpProposal.getDomain().imposeBound(y, Variable.TSTAR, new OmegaBigDecimal(BigDecimal.valueOf(50)));
        arMultivariateExpProposal.getDomain().imposeBound(z, Variable.TSTAR, new OmegaBigDecimal(BigDecimal.valueOf(50)));


        return arMultivariateExpProposal;
    }

    private static GEN getARMultivariateUniformProposal(Variable x, Variable y, Variable z) {
        //Multivariate uniform proposal construction for the AR method
        Expolynomial expolynomialUniformProposalFunction = new Expolynomial();

        Exmonomial exUniformProposal = new Exmonomial(BigDecimal.valueOf(0.2));
        exUniformProposal.addAtomicTerm(new MonomialTerm(x, 0));
        Exmonomial ex2UniformProposal = new Exmonomial(BigDecimal.valueOf(0.2));
        ex2UniformProposal.addAtomicTerm(new MonomialTerm(y, 0));
        Exmonomial ex3UniformProposal = new Exmonomial(BigDecimal.valueOf(0.1));
        ex3UniformProposal.addAtomicTerm(new MonomialTerm(z, 0));
        exUniformProposal.multiply(ex2UniformProposal);
        exUniformProposal.multiply(ex3UniformProposal);
        expolynomialUniformProposalFunction.addExmonomial(exUniformProposal);

        GEN arMultivariateUniformProposal = new GEN(new DBMZone(x,y,z), expolynomialUniformProposalFunction);
        arMultivariateUniformProposal.getDomain().imposeBound(x, Variable.TSTAR, new OmegaBigDecimal(BigDecimal.valueOf(5))); // upper bound
        arMultivariateUniformProposal.getDomain().imposeBound(Variable.TSTAR, x, new OmegaBigDecimal(BigDecimal.valueOf(0)).negate()); // lower bound
        arMultivariateUniformProposal.getDomain().imposeBound(y, Variable.TSTAR, new OmegaBigDecimal(BigDecimal.valueOf(5))); // upper bound
        arMultivariateUniformProposal.getDomain().imposeBound(Variable.TSTAR, y, new OmegaBigDecimal(BigDecimal.valueOf(0)).negate()); // lower bound
        arMultivariateUniformProposal.getDomain().imposeBound(z, Variable.TSTAR, new OmegaBigDecimal(BigDecimal.valueOf(10))); // upper bound
        arMultivariateUniformProposal.getDomain().imposeBound(Variable.TSTAR, z, new OmegaBigDecimal(BigDecimal.valueOf(0)).negate()); // lower bound
        return arMultivariateUniformProposal;
    }

}



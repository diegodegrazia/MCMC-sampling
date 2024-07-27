package org.example;

import org.example.Plots.ErlangPlot;
import org.example.Plots.TruncatedExpPlot;
import org.example.Plots.BivariatePlot;
import org.oristool.math.OmegaBigDecimal;
import org.oristool.math.expression.Exmonomial;
import org.oristool.math.expression.Expolynomial;
import org.oristool.math.expression.MonomialTerm;
import org.oristool.math.expression.Variable;
import org.oristool.math.function.EXP;
import org.oristool.math.function.Erlang;
import org.oristool.math.function.GEN;

import javax.swing.*;
import java.math.BigDecimal;
import java.math.MathContext;
import java.util.*;

public class Main {
    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            Scanner scanner = new Scanner(System.in);
            int k = 0;
            BigDecimal lambda = BigDecimal.ZERO;
            OmegaBigDecimal leftTruncatedExpBound;
            OmegaBigDecimal rightTruncatedExpBound;
            int numSamples = 0;
            // TODO: Ask the proposal ar function to the user and pass it as parameter
            GEN uniformProposal = GEN.newUniform(OmegaBigDecimal.ZERO, new OmegaBigDecimal(BigDecimal.valueOf(7)));
            double c = 1.31; // constant c for the ar method
            int[] arHitsAndMisses = new int[2]; // 0 for the hits, 1 for the misses
            int[] mhHitsAndMisses = new int[2]; // 0 for the hits, 1 for the misses

            Variable x = new Variable("x");
            Variable y = new Variable("y");
            /*Expolynomial myFunction = Expolynomial.fromString("4.0*x");
            System.out.println(myFunction.getVariables());
            MonomialTerm m = new MonomialTerm(y, 1);
            Exmonomial ex = new Exmonomial(BigDecimal.valueOf(3.0));
            ex.addAtomicTerm(m);
            myFunction.addExmonomial(ex);
            System.out.println(myFunction.getExmonomials()); */

            // Input shape parameter
            while (k <= 0) {
                System.out.print("Enter the shape parameter (k > 0): ");
                try {
                    k = scanner.nextInt();
                    if (k <= 0) {
                        System.out.println("Shape parameter must be greater than 0.");
                    }
                } catch (InputMismatchException e) {
                    System.out.println("Invalid input. Please enter a positive integer.");
                    scanner.next(); // Clear invalid input
                }
            }

            // Input rate parameter
            while (lambda.doubleValue() <= 0) {
                System.out.print("Enter the rate parameter (lambda > 0): ");
                try {
                    lambda = BigDecimal.valueOf(scanner.nextDouble());
                    if (lambda.doubleValue() <= 0) {
                        System.out.println("Rate parameter must be greater than 0.");
                    }
                } catch (InputMismatchException e) {
                    System.out.println("Invalid input. Please enter a positive number.");
                    scanner.next(); // Clear invalid input
                }
            }

            // Input leftTruncatedExpBound
            while (true) {
                System.out.print("Enter the support leftTruncatedExpBound (leftTruncatedExpBound >= 0): ");
                try {
                    BigDecimal eftInput = BigDecimal.valueOf(scanner.nextDouble());
                    if (eftInput.compareTo(BigDecimal.ZERO) < 0) {
                        System.out.println("leftTruncatedExpBound must be greater than or equal to 0.");
                    } else {
                        leftTruncatedExpBound = new OmegaBigDecimal(eftInput);
                        break;
                    }
                } catch (InputMismatchException e) {
                    System.out.println("Invalid input. Please enter a non-negative number.");
                    scanner.next(); // Clear invalid input
                }
            }

            // Input rightTruncatedExpBound
            while (true) {
                System.out.print("Enter the support rightTruncatedExpBound (rightTruncatedExpBound > leftTruncatedExpBound): ");
                try {
                    BigDecimal lftInput = BigDecimal.valueOf(scanner.nextDouble());
                    if (lftInput.compareTo(leftTruncatedExpBound.bigDecimalValue()) <= 0) {
                        System.out.println("rightTruncatedExpBound must be greater than leftTruncatedExpBound.");
                    } else {
                        rightTruncatedExpBound = new OmegaBigDecimal(lftInput);
                        break;
                    }
                } catch (InputMismatchException e) {
                    System.out.println("Invalid input. Please enter a number greater than leftTruncatedExpBound.");
                    scanner.next(); // Clear invalid input
                }
            }

            // Input number of samples
            while (numSamples <= 0) {
                System.out.print("Enter the number of samples for Metropolis-Hastings (positive integer): ");
                try {
                    numSamples = scanner.nextInt();
                    if (numSamples <= 0) {
                        System.out.println("Number of samples must be greater than 0.");
                    }
                } catch (InputMismatchException e) {
                    System.out.println("Invalid input. Please enter a positive integer.");
                    scanner.next(); // Clear invalid input
                }
            }

            // Create Erlang distribution
            Erlang erlang = new Erlang(x, k, lambda);
            Erlang erlang2 = new Erlang(y, k, lambda);
            GEN bivariateIID = erlang.cartesianProduct(erlang2);
            GEN boundedExp = GEN.newTruncatedExp(x, lambda, leftTruncatedExpBound, rightTruncatedExpBound); // a simple truncated exp to be plotted
            EXP exp = new EXP(x, lambda); // for the construction of the Erlang as sum of k exponentials


            // Metropolis-Hastings sampling
            double[][] bivariateSamples = metropolisHastingsMultivariate(bivariateIID, numSamples);
            double[] mhSamples = metropolisHastings(erlang, numSamples, mhHitsAndMisses);

            // Acceptance-Rejection sampling
            // need a list instead of an array because I don't want a fixed length,
            // the plot would be condensed to 0 values for high misses. So I add only hits in the list.
            List<Double> arSamples;
            arSamples = acceptanceRejection(erlang, numSamples, uniformProposal, c, arHitsAndMisses);

            // Sum of k Exponentials
            double[] sumOfExponentialsSamples = generateSumOfExponentialsSamples(exp, k, numSamples);

            // Plot the histograms and PDFs
            ErlangPlot.plotCharts(mhSamples, erlang, mhHitsAndMisses[0], mhHitsAndMisses[1]);
            TruncatedExpPlot.plotCharts(boundedExp); // No need for hits or misses beacuase I apply inv(CDF).
            BivariatePlot.plotBivariateCharts(bivariateSamples);
            ErlangPlot.plotAcceptanceRejectionCharts(arSamples, erlang, uniformProposal, c, arHitsAndMisses[0], arHitsAndMisses[1]);
            ErlangPlot.plotSumOfExponentialsCharts(sumOfExponentialsSamples, erlang);
        });
    }

    private static double[] metropolisHastings(Erlang targetDistribution, int numSamples, int[] mhHitsAndMisses) {
        double[] samples = new double[numSamples];
        Random random = new Random();
        MathContext mc = MathContext.DECIMAL128;
        OmegaBigDecimal current = new OmegaBigDecimal(new BigDecimal(random.nextDouble(), mc));
        OmegaBigDecimal proposal;
        OmegaBigDecimal alpha;

        Map<Variable, OmegaBigDecimal> variableMap = new HashMap<>();

        for (int i = 0; i < numSamples; i++) {
            proposal = current.add(new OmegaBigDecimal(new BigDecimal(random.nextGaussian(), mc)));

            if (proposal.compareTo(OmegaBigDecimal.ZERO) <= 0) { // if proposal <= 0 then I discard it and stay in the current state.
                samples[i] = current.doubleValue();
                mhHitsAndMisses[1]++;
                continue;
            }

            variableMap.put(targetDistribution.getVariable(), current);
            OmegaBigDecimal currentDensityValue = targetDistribution.getDensity().evaluate(variableMap);

            variableMap.put(targetDistribution.getVariable(), proposal);
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

    private static List<Double> acceptanceRejection(Erlang targetDistribution, int numSamples, GEN uniformProposal, double c, int [] hitsAndMisses) {
        List<Double> samples = new ArrayList<>();
        Random random = new Random();

        Map<Variable, OmegaBigDecimal> variableMap = new HashMap<>();
        Map<Variable, OmegaBigDecimal> uniformVariableMap = new HashMap<>();
        // TODO: Try to refactor this horrible way to get the variable
        Variable[] vars = uniformProposal.getDomain().getVariables().toArray(new Variable[uniformProposal.getDomain().getVariables().size()]);
        Variable proposaleVariable = vars[1]; // The uniform proposal is supposed to be a monovariate, vars[0] is the DBM reference variable x_0

        for (int i = 0; i < numSamples; i++) {

            double u = random.nextDouble();
            // uniform proposal sample y
            double y = uniformProposal.getDomainsEFT().doubleValue() + random.nextDouble() * (uniformProposal.getDomainsLFT().doubleValue() - uniformProposal.getDomainsEFT().doubleValue());
            variableMap.put(targetDistribution.getVariable(), new OmegaBigDecimal(BigDecimal.valueOf(y)));
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

    private static double[][] metropolisHastingsMultivariate(GEN multivariate, int numSamples) {

        Variable[] vars = multivariate.getDomain().getVariables().toArray(new Variable[multivariate.getDomain().getVariables().size() - 1]);
        List<Variable> variables = new ArrayList<>(Arrays.asList(vars));
        variables.remove(0);
        Random random = new Random();
        double[][] samples = new double[numSamples][variables.size()];
        MathContext mc = MathContext.DECIMAL128;
        OmegaBigDecimal[] currents = new OmegaBigDecimal[variables.size()];
        for (int i = 0; i < variables.size(); i++) {
            currents[i] = new OmegaBigDecimal(new BigDecimal(random.nextDouble(), mc));
        }
        OmegaBigDecimal[] proposals = new OmegaBigDecimal[variables.size()];
        OmegaBigDecimal alpha;

        Map<Variable, OmegaBigDecimal>[] variableMaps = new HashMap[2]; // for current and proposal
        variableMaps[0] = new HashMap<>(); // 0 for current
        variableMaps[1] = new HashMap<>(); // 1 for proposal

        outerLoop:
        for (int i = 0; i < numSamples; i++) {
            // Generate proposal by adding a small random perturbation
            for (int j = 0; j < variables.size(); j++) {
                proposals[j] = currents[j].add(new OmegaBigDecimal(new BigDecimal(random.nextGaussian(), mc)));
                // Ensure proposal is greater than zero
                if (proposals[j].compareTo(OmegaBigDecimal.ZERO) <= 0) {
                    for (int k = 0; k < variables.size(); k++) {
                        samples[i][j] = currents[j].doubleValue();
                    }
                    continue outerLoop;
                }
                // Preparing the variable map for current and proposal density evaluation
                variableMaps[0].put(variables.get(j), currents[j]);
                variableMaps[1].put(variables.get(j), proposals[j]);

            }
            OmegaBigDecimal currentDensityValue = multivariate.getDensity().evaluate(variableMaps[0]);

            OmegaBigDecimal proposalDensityValue = multivariate.getDensity().evaluate(variableMaps[1]);

            // Calculate acceptance ratio alpha
            alpha = proposalDensityValue.divide(currentDensityValue.bigDecimalValue(), mc).min(OmegaBigDecimal.ONE);

            // Accept or reject proposal
            if (new BigDecimal(random.nextDouble(), mc).compareTo(alpha.bigDecimalValue()) < 0) {
                for (int j = 0; j < variables.size(); j++)
                    currents[j] = proposals[j];
            }
            for (int j = 0; j < variables.size(); j++)
                samples[i][j] = currents[j].doubleValue();
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
}

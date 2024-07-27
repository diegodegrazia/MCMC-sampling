package org.example.Plots;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RefineryUtilities;
import org.oristool.math.OmegaBigDecimal;
import org.oristool.math.expression.Variable;
import org.oristool.math.function.Erlang;
import org.oristool.math.function.GEN;

import javax.swing.*;
import java.awt.*;
import java.math.BigDecimal;
import java.math.MathContext;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ErlangPlot extends ApplicationFrame {

    public ErlangPlot(String title, double[] samples, Erlang erlang, boolean includeProposal, int nHits, int nMisses) { // For the MH plot
        super(title);
        JFreeChart histogram = createHistogramChart(createDataset(samples), nHits, nMisses);
        JFreeChart pdfChart = createPDFChart(erlang, includeProposal, null, 0);

        ChartPanel histogramPanel = new ChartPanel(histogram);
        histogramPanel.setPreferredSize(new java.awt.Dimension(800, 400));
        ChartPanel pdfPanel = new ChartPanel(pdfChart);
        pdfPanel.setPreferredSize(new java.awt.Dimension(800, 400));

        JPanel panel = new JPanel(new GridLayout(2, 1));
        panel.add(histogramPanel);
        panel.add(pdfPanel);

        setContentPane(panel);
    }

    public ErlangPlot(String title, double[] samples, Erlang erlang, boolean includeProposal, GEN uniformProposal, double c, int nHits, int nMisses) { // For the AR plot
        super(title);

        JFreeChart histogram = createHistogramChart(createDataset(samples), nHits, nMisses);
        JFreeChart pdfChart = createPDFChart(erlang, includeProposal, uniformProposal, c);

        ChartPanel histogramPanel = new ChartPanel(histogram);
        histogramPanel.setPreferredSize(new java.awt.Dimension(800, 400));
        ChartPanel pdfPanel = new ChartPanel(pdfChart);
        pdfPanel.setPreferredSize(new java.awt.Dimension(800, 400));

        JPanel panel = new JPanel(new GridLayout(2, 1));
        panel.add(histogramPanel);
        panel.add(pdfPanel);

        setContentPane(panel);
    }

    private HistogramDataset createDataset(double[] samples) {
        HistogramDataset dataset = new HistogramDataset();
        dataset.addSeries("Samples", samples, 80); // 80 bins
        return dataset;
    }

    private JFreeChart createHistogramChart(HistogramDataset dataset, int nHits, int nMisses) {
        JFreeChart histogramChart = ChartFactory.createHistogram(
                "Histogram of Samples",
                "Value",
                "Frequency",
                dataset,
                PlotOrientation.VERTICAL,
                true,
                true,
                false
        );

        XYPlot plot = (XYPlot) histogramChart.getPlot();
        NumberAxis xAxis = (NumberAxis) plot.getDomainAxis();
        xAxis.setRange(0, 50);

        // Add hits, misses, and efficiency to the chart title
        String chartTitle = String.format("Histogram of Samples (Hits: %d, Misses: %d, Efficiency: %.2f%%)",nHits, nMisses, ((float)nHits/ (nHits+(float)nMisses))*100);
        histogramChart.setTitle(chartTitle);

        return histogramChart;
    }

    private JFreeChart createPDFChart(Erlang erlang, boolean includeProposal, GEN uniformProposal, double c) {
        JFreeChart pdfChart = ChartFactory.createXYLineChart(
                "Erlang PDF",
                "x",
                "Probability Density",
                createPDFDataset(erlang, includeProposal, uniformProposal, c),
                PlotOrientation.VERTICAL,
                true,
                true,
                false
        );

        XYPlot plot = (XYPlot) pdfChart.getPlot();
        NumberAxis xAxis = (NumberAxis) plot.getDomainAxis();
        xAxis.setRange(0, 50);

        return pdfChart;
    }

    private XYSeriesCollection createPDFDataset(Erlang erlang, boolean includeProposal, GEN uniformProposal, double c) {
        XYSeries erlangSeries = new XYSeries("Erlang PDF");
        XYSeries proposalSeries = new XYSeries("Proposal Distribution");

        double start = 0;
        double end = 50;
        int numPoints = 10000;
        double step = (end - start) / numPoints;

        for (double x = start; x <= end; x += step) {
            double pdfValue = calculateErlangPDF(erlang, x);
            erlangSeries.add(x, pdfValue);
        }

        if (includeProposal) {
            for (double x = uniformProposal.getDomainsEFT().doubleValue(); x <= uniformProposal.getDomainsLFT().doubleValue(); x += step){
                double proposalValue = calculateProposalPDF(x, uniformProposal, c);
                proposalSeries.add(x, proposalValue);
            }
        }

        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(erlangSeries);
        if (includeProposal) {
            dataset.addSeries(proposalSeries);
        }

        return dataset;
    }

    private double calculateErlangPDF(Erlang erlang, double x) {
        Map<Variable, OmegaBigDecimal> variableMap = new HashMap<>();
        variableMap.put(erlang.getVariable(), new OmegaBigDecimal(new BigDecimal(x, MathContext.DECIMAL128)));
        OmegaBigDecimal pdfValue = erlang.getDensity().evaluate(variableMap);
        return pdfValue.doubleValue();
    }

    private double calculateProposalPDF(double x, GEN uniformProposal, double c) {
        if (uniformProposal == null) {
            return 0.0;
        }

        Map<Variable, OmegaBigDecimal> variableMap = new HashMap<>();
        Variable[] vars = uniformProposal.getDomain().getVariables().toArray(new Variable[0]);
        Variable proposalVariable = vars[1]; // Assuming uniform proposal is univariate
        variableMap.put(proposalVariable, new OmegaBigDecimal(new BigDecimal(x, MathContext.DECIMAL128)));

        OmegaBigDecimal proposalDensityValue = uniformProposal.getDensity().evaluate(variableMap);
        return proposalDensityValue.doubleValue() * c;
    }

    public static void plotCharts(double[] samples, Erlang erlang, int nHits, int nMisses) {
        ErlangPlot chart = new ErlangPlot("Metropolis-Hastings and Erlang PDF", samples, erlang, false, nHits, nMisses);
        chart.pack();
        RefineryUtilities.centerFrameOnScreen(chart);
        chart.setVisible(true);
    }

    public static void plotAcceptanceRejectionCharts(List<Double> samples, Erlang erlang, GEN uniformProposal, double c, int nHits, int nMisses) {
        double[] arSample = new double[samples.size()];
        for (int i = 0; i < samples.size(); i++) {
            arSample[i] = samples.get(i);
        }
        ErlangPlot chart = new ErlangPlot("Acceptance-Rejection and Erlang PDF", arSample, erlang, true, uniformProposal, c, nHits, nMisses);
        chart.pack();
        RefineryUtilities.centerFrameOnScreen(chart);
        chart.setVisible(true);
    }

    public static void plotSumOfExponentialsCharts(double[] samples, Erlang erlang) {
        // TODO: I don't really like these useless 0s. Check again the way you eluded the lack of default parameters in java
        ErlangPlot chart = new ErlangPlot("Sum of Exponentials and Erlang PDF", samples, erlang, false, 0, 0);
        chart.pack();
        RefineryUtilities.centerFrameOnScreen(chart);
        chart.setVisible(true);
    }
}

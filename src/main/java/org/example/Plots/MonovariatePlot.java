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
import org.oristool.math.function.GEN;
import org.oristool.math.function.PartitionedGEN;

import javax.swing.*;
import java.awt.*;
import java.math.BigDecimal;
import java.math.MathContext;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class MonovariatePlot extends ApplicationFrame {

    public MonovariatePlot(String title, double[] samples, GEN genFunction, PartitionedGEN piecewiseProposal, double c, int nHits, int nMisses) {
        super(title);

        JFreeChart histogram = createHistogramChart(createDataset(samples), nHits, nMisses);
        JFreeChart pdfChart = createPiecewisePDFChart(genFunction, piecewiseProposal, c);

        ChartPanel histogramPanel = new ChartPanel(histogram);
        histogramPanel.setPreferredSize(new java.awt.Dimension(800, 400));
        ChartPanel pdfPanel = new ChartPanel(pdfChart);
        pdfPanel.setPreferredSize(new java.awt.Dimension(800, 400));

        JPanel panel = new JPanel(new GridLayout(2, 1));
        panel.add(histogramPanel);
        panel.add(pdfPanel);

        setContentPane(panel);
    }

    public MonovariatePlot(String title, double[] samples, GEN genFunction, boolean includeProposal, int nHits, int nMisses) { // For the MH plot
        super(title);
        JFreeChart histogram = createHistogramChart(createDataset(samples), nHits, nMisses);
        JFreeChart pdfChart = createPDFChart(genFunction, includeProposal, null, 0);

        ChartPanel histogramPanel = new ChartPanel(histogram);
        histogramPanel.setPreferredSize(new java.awt.Dimension(800, 400));
        ChartPanel pdfPanel = new ChartPanel(pdfChart);
        pdfPanel.setPreferredSize(new java.awt.Dimension(800, 400));

        JPanel panel = new JPanel(new GridLayout(2, 1));
        panel.add(histogramPanel);
        panel.add(pdfPanel);

        setContentPane(panel);
    }

    public MonovariatePlot(String title, double[] samples, GEN genFunction, boolean includeProposal, GEN uniformProposal, double c, int nHits, int nMisses) { // For the AR plot
        super(title);

        JFreeChart histogram = createHistogramChart(createDataset(samples), nHits, nMisses);
        JFreeChart pdfChart = createPDFChart(genFunction, includeProposal, uniformProposal, c);

        ChartPanel histogramPanel = new ChartPanel(histogram);
        histogramPanel.setPreferredSize(new java.awt.Dimension(800, 400));
        ChartPanel pdfPanel = new ChartPanel(pdfChart);
        pdfPanel.setPreferredSize(new java.awt.Dimension(800, 400));

        JPanel panel = new JPanel(new GridLayout(2, 1));
        panel.add(histogramPanel);
        panel.add(pdfPanel);

        setContentPane(panel);
    }

    public MonovariatePlot(String title, GEN genFunction) {
        super(title);
        JFreeChart pdfChart = createPDFChart(genFunction, false, null, 0);

        ChartPanel pdfPanel = new ChartPanel(pdfChart);
        pdfPanel.setPreferredSize(new java.awt.Dimension(800, 400));

        JPanel panel = new JPanel(new BorderLayout());
        panel.add(pdfPanel, BorderLayout.CENTER);

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
        String chartTitle = String.format("Histogram of Samples (Hits: %d, Misses: %d, Efficiency: %.2f%%)", nHits, nMisses, ((float) nHits / (nHits + (float) nMisses)) * 100);
        histogramChart.setTitle(chartTitle);

        return histogramChart;
    }

    private JFreeChart createPDFChart(GEN genFunction, boolean includeProposal, GEN uniformProposal, double c) {
        JFreeChart pdfChart = ChartFactory.createXYLineChart(
                "PDF",
                "x",
                "Probability Density",
                createPDFDataset(genFunction, includeProposal, uniformProposal, c),
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

    private XYSeriesCollection createPDFDataset(GEN genFunction, boolean includeProposal, GEN uniformProposal, double c) {
        XYSeries pdfSeries = new XYSeries("PDF");
        XYSeries proposalSeries = new XYSeries("Proposal Distribution");

        double start = genFunction.getDomainsEFT().doubleValue();
        double end = genFunction.getDomainsLFT().doubleValue();
        if (end == OmegaBigDecimal.POSITIVE_INFINITY.doubleValue())
            end = 50;
        int numPoints = 10000;
        double step = (end - start) / numPoints;

        for (double x = start; x <= end; x += step) {
            double pdfValue = calculatePDF(genFunction, x);
            pdfSeries.add(x, pdfValue);
        }

        if (includeProposal) {
            for (double x = uniformProposal.getDomainsEFT().doubleValue(); x <= uniformProposal.getDomainsLFT().doubleValue(); x += step) {
                double proposalValue = calculateProposalPDF(x, uniformProposal, c);
                proposalSeries.add(x, proposalValue);
            }
        }

        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(pdfSeries);
        if (includeProposal) {
            dataset.addSeries(proposalSeries);
        }

        return dataset;
    }

    private JFreeChart createPiecewisePDFChart(GEN genFunction, PartitionedGEN piecewiseProposal, double c) {
        JFreeChart pdfChart = ChartFactory.createXYLineChart(
                "PDF with Piecewise Proposal",
                "x",
                "Probability Density",
                createPiecewisePDFDataset(genFunction, piecewiseProposal, c),
                PlotOrientation.VERTICAL,
                true,
                true,
                false
        );

        XYPlot plot = (XYPlot) pdfChart.getPlot();
        NumberAxis xAxis = (NumberAxis) plot.getDomainAxis();
        xAxis.setRange(0, 50);

        // Set a single color for the combined piecewise series
        plot.getRenderer().setSeriesPaint(1, Color.BLUE); // Assuming the piecewise proposal is the second series

        return pdfChart;
    }

    private XYSeriesCollection createPiecewisePDFDataset(GEN genFunction, PartitionedGEN piecewiseProposal, double c) {
        XYSeries pdfSeries = new XYSeries("PDF");
        XYSeries proposalSeries = new XYSeries("Piecewise Proposal");

        double start = genFunction.getDomainsEFT().doubleValue();
        double end = genFunction.getDomainsLFT().doubleValue();
        if (end == OmegaBigDecimal.POSITIVE_INFINITY.doubleValue()) end = 50;
        int numPoints = 10000;
        double step = (end - start) / numPoints;

        // Create PDF series
        for (double x = start; x <= end; x += step) {
            double pdfValue = calculatePDF(genFunction, x);
            pdfSeries.add(x, pdfValue);
        }

        // Combine proposal series into one
        if (piecewiseProposal != null) {
            for (GEN piece : piecewiseProposal.getFunctions()) {
                double pieceStart = piece.getDomainsEFT().doubleValue();
                double pieceEnd = piece.getDomainsLFT().doubleValue();
                if (pieceEnd == OmegaBigDecimal.POSITIVE_INFINITY.doubleValue()) pieceEnd = 50;

                for (double x = pieceStart; x <= pieceEnd; x += step) {
                    double proposalValue = calculateProposalPDF(x, piece, c);
                    proposalSeries.add(x, proposalValue);
                }
            }
        }

        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(pdfSeries);
        if (piecewiseProposal != null) {
            dataset.addSeries(proposalSeries); // Add the combined series to the dataset
        }

        return dataset;
    }

    private double calculatePDF(GEN genFunction, double x) {
        Variable[] vars = genFunction.getDomain().getVariables().toArray(new Variable[0]);
        Variable variable = vars[1];
        Map<Variable, OmegaBigDecimal> variableMap = new HashMap<>();
        variableMap.put(variable, new OmegaBigDecimal(new BigDecimal(x, MathContext.DECIMAL128)));
        OmegaBigDecimal pdfValue = genFunction.getDensity().evaluate(variableMap);
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


    public static void plotCharts(double[] samples, GEN genFunction, int nHits, int nMisses) {
        MonovariatePlot chart = new MonovariatePlot("Metropolis-Hastings and PDF", samples, genFunction, false, nHits, nMisses);
        chart.pack();
        RefineryUtilities.centerFrameOnScreen(chart);
        chart.setVisible(true);
    }

    public static void plotCharts(GEN genFunction) {
        MonovariatePlot chart = new MonovariatePlot("Metropolis-Hastings and PDF", genFunction);
        chart.pack();
        RefineryUtilities.centerFrameOnScreen(chart);
        chart.setVisible(true);
    }

    public static void plotAcceptanceRejectionCharts(List<Double> samples, GEN genFunction, GEN uniformProposal, double c, int nHits, int nMisses) {
        double[] arSample = new double[samples.size()];
        for (int i = 0; i < samples.size(); i++) {
            arSample[i] = samples.get(i);
        }
        MonovariatePlot chart = new MonovariatePlot("Acceptance-Rejection and PDF", arSample, genFunction, true, uniformProposal, c, nHits, nMisses);
        chart.pack();
        RefineryUtilities.centerFrameOnScreen(chart);
        chart.setVisible(true);
    }

    public static void plotSumOfExponentialsCharts(double[] samples, GEN genFunction) {
        MonovariatePlot chart = new MonovariatePlot("Sum of Exponentials and PDF", samples, genFunction, false, 0, 0);
        chart.pack();
        RefineryUtilities.centerFrameOnScreen(chart);
        chart.setVisible(true);
    }

    public static void plotPiecewiseAcceptanceRejectionCharts(List<Double> samples, GEN genFunction, PartitionedGEN piecewiseProposal, double c, int nHits, int nMisses) {
        double[] arSample = new double[samples.size()];
        for (int i = 0; i < samples.size(); i++) {
            arSample[i] = samples.get(i);
        }
        MonovariatePlot chart = new MonovariatePlot("Acceptance-Rejection with Piecewise Proposal", arSample, genFunction, piecewiseProposal, c, nHits, nMisses);
        chart.pack();
        RefineryUtilities.centerFrameOnScreen(chart);
        chart.setVisible(true);
    }
}

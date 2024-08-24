package org.example.Plots;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.xy.DefaultXYDataset;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;
import java.util.List;

public class MultivariatePlot {

    public static void plotMultivariateCharts(List<double[]> samplesList, int nHits, int nMisses, String algorithmName) {
        if (samplesList.isEmpty()) return;

        int dim = samplesList.get(0).length;

        List<JPanel> scatterPlots = new ArrayList<>();
        List<JPanel> histogramPlots = new ArrayList<>();

        // Prepare histograms
        for (int i = 0; i < dim; i++) {
            histogramPlots.add(createHistogram(samplesList, i, nHits, nMisses, algorithmName));
        }

        // Prepare scatter plots
        for (int i = 0; i < dim - 1; i++) {
            for (int j = i + 1; j < dim; j++) {
                scatterPlots.add(createScatterPlot(samplesList, i, j, nHits, nMisses, algorithmName));
            }
        }

        displayChartsIndividually(scatterPlots, histogramPlots, algorithmName);
    }

    private static JPanel createScatterPlot(List<double[]> samplesList, int xIndex, int yIndex, int nHits, int nMisses, String algorithmName) {
        DefaultXYDataset dataset = new DefaultXYDataset();

        double[][] data = new double[2][samplesList.size()];
        for (int i = 0; i < samplesList.size(); i++) {
            data[0][i] = samplesList.get(i)[xIndex];
            data[1][i] = samplesList.get(i)[yIndex];
        }

        dataset.addSeries("Scatter Plot", data);

        String title = String.format("%s: Scatter Plot (Variable %d vs %d) (nHits: %d, nMisses: %d, Efficiency: %.2f%%)",
                algorithmName, xIndex + 1, yIndex + 1, nHits, nMisses, (double) nHits / (nHits + nMisses) * 100);

        JFreeChart scatterPlot = ChartFactory.createScatterPlot(
                title,
                "Variable " + (xIndex + 1),
                "Variable " + (yIndex + 1),
                dataset,
                PlotOrientation.VERTICAL,
                true,
                true,
                false
        );

        XYPlot plot = scatterPlot.getXYPlot();
        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer(false, true);
        plot.setRenderer(renderer);

        NumberAxis domain = (NumberAxis) plot.getDomainAxis();
        domain.setRange(domain.getRange());

        NumberAxis range = (NumberAxis) plot.getRangeAxis();
        range.setRange(range.getRange());

        return new ChartPanel(scatterPlot);
    }

    private static JPanel createHistogram(List<double[]> samplesList, int varIndex, int nHits, int nMisses, String algorithmName) {
        HistogramDataset dataset = new HistogramDataset();
        double[] data = new double[samplesList.size()];
        for (int i = 0; i < samplesList.size(); i++) {
            data[i] = samplesList.get(i)[varIndex];
        }

        dataset.addSeries("Variable " + (varIndex + 1), data, 50);

        String title = String.format("%s: Histogram (Variable %d) (nHits: %d, nMisses: %d, Efficiency: %.2f%%)",
                algorithmName, varIndex + 1, nHits, nMisses, (double) nHits / (nHits + nMisses) * 100);

        JFreeChart histogram = ChartFactory.createHistogram(
                title,
                "Value",
                "Frequency",
                dataset,
                PlotOrientation.VERTICAL,
                true,
                true,
                false
        );

        return new ChartPanel(histogram);
    }

    private static void displayChartsIndividually(List<JPanel> scatterPlots, List<JPanel> histogramPlots, String algorithmName) {
        // Display each histogram in its own window
        for (int i = 0; i < histogramPlots.size(); i++) {
            JFrame frame = new JFrame(algorithmName + " - Histogram for Variable " + (i + 1));
            frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
            frame.add(histogramPlots.get(i), BorderLayout.CENTER);
            frame.pack();
            frame.setVisible(true);
        }

        // Display each scatter plot in its own window
        int count = 1;
        for (JPanel chart : scatterPlots) {
            JFrame frame = new JFrame(algorithmName + " - Scatter Plot " + count++);
            frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
            frame.add(chart, BorderLayout.CENTER);
            frame.pack();
            frame.setVisible(true);
        }
    }
}

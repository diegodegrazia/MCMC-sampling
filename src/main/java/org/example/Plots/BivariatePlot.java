package org.example.Plots;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RefineryUtilities;

import javax.swing.*;
import java.awt.*;

public class BivariatePlot extends ApplicationFrame {

    public BivariatePlot(String title, double[][] samples) {
        super(title);
        JFreeChart scatterPlot = createScatterPlot(samples);

        ChartPanel chartPanel = new ChartPanel(scatterPlot);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));

        setContentPane(chartPanel);
    }

    private JFreeChart createScatterPlot(double[][] samples) {
        XYSeries series = new XYSeries("Bivariate Samples");

        for (double[] sample : samples) {
            series.add(sample[0], sample[1]);
        }

        XYSeriesCollection dataset = new XYSeriesCollection(series);

        JFreeChart chart = ChartFactory.createScatterPlot(
                "Bivariate Scatter Plot",
                "X",
                "Y",
                dataset,
                PlotOrientation.VERTICAL,
                true,
                true,
                false
        );

        XYPlot plot = (XYPlot) chart.getPlot();
        NumberAxis domainAxis = (NumberAxis) plot.getDomainAxis();
        NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
        domainAxis.setRange(0.0, 50.0); // Adjust the range as needed
        rangeAxis.setRange(0.0, 50.0);  // Adjust the range as needed

        return chart;
    }

    public static void plotBivariateCharts(double[][] samples) {
        BivariatePlot chart = new BivariatePlot("Bivariate Metropolis-Hastings Samples", samples);
        chart.pack();
        RefineryUtilities.centerFrameOnScreen(chart);
        chart.setVisible(true);
    }
}

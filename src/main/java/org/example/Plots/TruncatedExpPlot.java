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
import org.oristool.math.OmegaBigDecimal;
import org.oristool.math.expression.Variable;
import org.oristool.math.function.GEN;

import javax.swing.*;
import java.awt.*;
import java.math.BigDecimal;
import java.math.MathContext;
import java.util.Map;
import java.util.HashMap;

public class TruncatedExpPlot extends ApplicationFrame {

    public TruncatedExpPlot(String title, GEN boundedExp) {
        super(title);
        JFreeChart pdfChart = createPDFChart(boundedExp);

        ChartPanel pdfPanel = new ChartPanel(pdfChart);
        pdfPanel.setPreferredSize(new java.awt.Dimension(800, 400));

        JPanel panel = new JPanel(new BorderLayout());
        panel.add(pdfPanel, BorderLayout.CENTER);

        setContentPane(panel);
    }

    private JFreeChart createPDFChart(GEN boundedExp) {
        JFreeChart pdfChart = ChartFactory.createXYLineChart(
                "Truncated Exponential PDF",
                "x",
                "Probability Density",
                createPDFDataset(boundedExp),
                PlotOrientation.VERTICAL,
                true,
                true,
                false
        );

        // Set axis range
        XYPlot plot = (XYPlot) pdfChart.getPlot();
        NumberAxis xAxis = (NumberAxis) plot.getDomainAxis();
        xAxis.setRange(boundedExp.getDomainsEFT().doubleValue(), boundedExp.getDomainsLFT().doubleValue());

        return pdfChart;
    }

    private XYSeriesCollection createPDFDataset(GEN boundedExp) {
        XYSeries series = new XYSeries("Truncated Exponential PDF");

        // Define the range for x values
        double start = boundedExp.getDomainsEFT().doubleValue();
        double end = boundedExp.getDomainsLFT().doubleValue();
        int numPoints = 1000; // Number of points for smoothness
        double step = (end - start) / numPoints;

        for (double x = start; x <= end; x += step) {
            double pdfValue = calculateExpPDF(boundedExp, x);
            series.add(x, pdfValue);
        }

        return new XYSeriesCollection(series);
    }

    private double calculateExpPDF(GEN boundedExp, double x) {
        Variable xVar = new Variable("x");
        Map<Variable, OmegaBigDecimal> variableMap = new HashMap<>();
        variableMap.put(xVar, new OmegaBigDecimal(new BigDecimal(x, MathContext.DECIMAL128)));
        OmegaBigDecimal pdfValue = boundedExp.getDensity().evaluate(variableMap);
        return pdfValue.doubleValue();
    }

    public static void plotCharts(GEN boundedExp) {
        TruncatedExpPlot chart = new TruncatedExpPlot("Truncated Exponential PDF", boundedExp);
        chart.pack();
        RefineryUtilities.centerFrameOnScreen(chart);
        chart.setVisible(true);
    }
}

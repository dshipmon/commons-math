package org.apache.commons.math4.analysis.integration.adaptive;

import org.apache.commons.math4.analysis.UnivariateFunction;
import org.apache.commons.math4.analysis.integration.UnivariateIntegrator;
import org.apache.commons.math4.analysis.integration.adaptive.node.AdaptiveQuadratureNode;
import org.apache.commons.math4.analysis.integration.gauss.GaussIntegrator;

import java.util.PriorityQueue;


public class AdaptiveIntegrator {
    PriorityQueue<AdaptiveQuadratureNode> integrationHeap = new PriorityQueue<>();
    private UnivariateFunction function;
    private UnivariateIntegrator integrator;

    public AdaptiveIntegrator(UnivariateFunction function, UnivariateIntegrator integrator){
        this.function = function;
        this.integrator = integrator;
    }

    public double integrate(double min, double max, int maxEval){
        Double Q = integrator.integrate(maxEval, function, min, max);
        return 0;
    }

    private double simpsonError(int a, int b){
        return 0;
    }
}

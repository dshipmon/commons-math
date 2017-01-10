package org.apache.commons.math4.analysis.integration.adaptive;

import org.apache.commons.math4.analysis.UnivariateFunction;
import org.apache.commons.math4.analysis.integration.UnivariateIntegrator;
import org.apache.commons.math4.analysis.integration.adaptive.node.AdaptiveQuadratureNode;
import org.apache.commons.math4.analysis.integration.gauss.GaussIntegrator;
import org.apache.commons.math4.analysis.integration.gauss.GaussIntegratorFactory;

import java.util.PriorityQueue;


public class AdaptiveIntegrator {
    PriorityQueue<AdaptiveQuadratureNode> integrationHeap;
    private UnivariateFunction function;

    public AdaptiveIntegrator(UnivariateFunction function, UnivariateIntegrator integrator){
        this.function = function;
    }

    public double integrate(double a, double b, double errorTolerance, int ruleOrder , int maxEval){
        GaussIntegratorFactory factory = new GaussIntegratorFactory();
        GaussIntegrator integratorLegendre = factory.legendre(ruleOrder, a, b);
        GaussIntegrator integratorKronrod = factory.kronrod(ruleOrder, a, b);
        double globalIntegral = integratorLegendre.integrate(function);
        double globalError = Math.abs(globalIntegral - integratorKronrod.integrate(function));

        if (globalError <= errorTolerance){
            return globalIntegral;
        }

        integrationHeap = new PriorityQueue<>(1);
        integrationHeap.add(new AdaptiveQuadratureNode(a, b, globalIntegral, globalError));

        int count = 0;
        double midpoint;
        double integralLeft;
        double integralRight;
        double errorLeft;
        double errorRight;
        AdaptiveQuadratureNode currentNode;
        while (globalError > errorTolerance && count < maxEval ){
            currentNode = integrationHeap.poll();
            midpoint = (currentNode.getA() + currentNode.getB()) / 2;

            integralLeft = (factory.legendre(ruleOrder, currentNode.getA(), midpoint)).integrate(function);
            integralRight = (factory.legendre(ruleOrder, midpoint, currentNode.getB())).integrate(function);

            errorLeft = Math.abs(integralLeft - (factory.kronrod(ruleOrder, currentNode.getA(), midpoint)).integrate(function));
            errorRight = Math.abs(integralRight - (factory.kronrod(ruleOrder, midpoint, currentNode.getB())).integrate(function));

            globalIntegral = globalIntegral - currentNode.getIntegral() + integralLeft + integralRight;
            globalError = globalError - currentNode.getError() + errorLeft + errorRight;

            integrationHeap.add(new AdaptiveQuadratureNode(currentNode.getA(), midpoint, integralLeft, errorLeft));
            integrationHeap.add(new AdaptiveQuadratureNode(midpoint, currentNode.getB(), integralRight, errorRight));
            count++;
        }
        return globalIntegral;
    }

    public double integrate(double a, double b, double errorTolerance){

        return 0;
    }

}

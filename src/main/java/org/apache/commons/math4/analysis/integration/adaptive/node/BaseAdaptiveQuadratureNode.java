package org.apache.commons.math4.analysis.integration.adaptive.node;


public abstract class BaseAdaptiveQuadratureNode<T> {
    private T a;
    private T b;
    private T integral;
    private T error;

    public BaseAdaptiveQuadratureNode(){}

    public BaseAdaptiveQuadratureNode(T a, T b, T integral, T error) {
        this.a = a;
        this.b = b;
        this.integral = integral;
        this.error = error;
    }

    public T getA() {
        return a;
    }

    public void setA(T a) {
        this.a = a;
    }

    public T getB() {
        return b;
    }

    public void setB(T b) {
        this.b = b;
    }

    public T getIntegral() {
        return integral;
    }

    public void setIntegral(T integral) {
        this.integral = integral;
    }

    public T getError() {
        return error;
    }

    public void setError(T error) {
        this.error = error;
    }
}

package org.apache.commons.math4.analysis.integration.adaptive.node;


public class AdaptiveQuadratureNode extends BaseAdaptiveQuadratureNode<Double> implements Comparable<AdaptiveQuadratureNode> {

    public AdaptiveQuadratureNode(Double a, Double b, Double approximateIntegral, Double error) {
        super(a, b, approximateIntegral, error);
    }

    @Override
    public int compareTo(AdaptiveQuadratureNode o) {
        return this.getError().compareTo(o.getError());
    }
}

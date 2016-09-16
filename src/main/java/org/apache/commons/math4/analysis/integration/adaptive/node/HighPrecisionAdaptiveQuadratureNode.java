package org.apache.commons.math4.analysis.integration.adaptive.node;

import java.math.BigDecimal;


public class HighPrecisionAdaptiveQuadratureNode extends BaseAdaptiveQuadratureNode<BigDecimal> implements Comparable<HighPrecisionAdaptiveQuadratureNode> {
    public HighPrecisionAdaptiveQuadratureNode() {
    }

    public HighPrecisionAdaptiveQuadratureNode(BigDecimal a, BigDecimal b, BigDecimal approximateIntegral, BigDecimal error) {
        super(a, b, approximateIntegral, error);
    }

    @Override
    public int compareTo(HighPrecisionAdaptiveQuadratureNode o) {
        return this.getError().compareTo(o.getError());
    }
}

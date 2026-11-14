## Start with two independent plots {#panels}

Use ggchord, dplyr and the original PAF in the [input bundle](../assets/tutorial-inputs.zip). Print course_center and course_detail independently, then connect one existing relationship.

```r
library(ggchord)
library(dplyr)
```

{{prepare:fusiondata}}

{{code:coursefusion}}

![Independent center plot](../assets/coursefusioncenter.png)

![Independent detail plot](../assets/coursefusiondetail.png)

![Two panels and one source relationship](../assets/coursefusion.png)

## Check intervals, panels and instances {#intervals}

Fusion is a composition, not a seventh Coord. Each panel retains its data, layers, scales and physical ratio. A side places a panel without rotating it. Endpoints must share units, be fully visible and route unambiguously; direct bp–aa links are unsupported.

```r
course_export <- export_ggchord_fusion(p_coursefusion, width = 12, height = 8)
course_export$connections
course_export$fronts
```

## Save at the same size {#save}

```r
ggsave("fusion-first.png", p_coursefusion,
  width = 12, height = 8, units = "in", dpi = 144)
```

Example: choose another real HSP and identify qpanel/spanel and its source row. Dotplot, facets and unregistered external plots are unsupported.

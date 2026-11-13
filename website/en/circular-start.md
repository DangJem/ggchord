## Length and a prepared annotation {#input}

Only lengths and prepared tables are needed. The ori row is the verified pUC19 annotation 867–1455, matching the original inspectorigin output; no DNA recognition is run. Dependencies: ggchord and dplyr.

```r
library(ggchord)
```

{{code:coursecircleinput}}

![Circular backbone](../assets/coursecircle.png)

## Add the feature and its label {#features-and-labels}

{{code:coursecirclefeature}}

![Circular feature and label](../assets/coursecirclefeature.png)

Locate both ori endpoints on the ticks and trace the label to anno. A circular display is not evidence of molecular topology.

## Exercise and save {#exercise}

Remove the label layer and rerun; backbone and source interval stay unchanged. Replace the row with your own known annotation.

```r
ggsave("circular-first.png", p_coursecirclefeature,
  width = 10, height = 8, units = "in", dpi = 144)
```

## 长度与已有注释 {#input}

只需要长度与已知表格。本课 ori 行来自已核对 pUC19 注释 867–1455，与原 inspectorigin 输出一致；不运行 DNA 识别。需要 ggchord、dplyr。

```r
library(ggchord)
```

{{code:coursecircleinput}}

![Circular backbone](../assets/coursecircle.png)

## 添加特征和标签 {#features-and-labels}

{{code:coursecirclefeature}}

![Circular feature and label](../assets/coursecirclefeature.png)

找到 ori 在刻度上的两端，确认文字来自 anno；圆形外观不表示已验证分子拓扑。

## 操作示例与保存 {#exercise}

去掉标签层再运行，骨架和源区间保持。用自己的已知注释替换该行。

```r
ggsave("circular-first.png", p_coursecirclefeature,
  width = 10, height = 8, units = "in", dpi = 144)
```

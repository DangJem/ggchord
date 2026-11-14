## 先独立画好两张图 {#panels}

需要 ggchord、dplyr 和[输入包](../assets/tutorial-inputs.zip)中的原始 PAF。先分别打印 course_center 和 course_detail，再组合一条已有关系。

```r
library(ggchord)
library(dplyr)
```

{{prepare:fusiondata}}

{{code:coursefusion}}

![独立中心图](../assets/coursefusioncenter.png)

![独立细节图](../assets/coursefusiondetail.png)

![Two panels and one source relationship](../assets/coursefusion.png)

## 核对区间、面板和实例 {#intervals}

Fusion 是组合对象，不是第七 Coord。每个面板保留自己的数据、图层、标度和比例；side 仅摆放，不旋转子图。两端须同单位、完整可见且路由唯一，不能直接连接 bp 与 aa。

```r
course_export <- export_ggchord_fusion(p_coursefusion, width = 12, height = 8)
course_export$connections
course_export$fronts
```

## 按相同尺寸保存 {#save}

```r
ggsave("fusion-first.png", p_coursefusion,
  width = 12, height = 8, units = "in", dpi = 144)
```

操作示例：换成另一条真实 HSP，指出 qpanel/spanel 与源行。Dotplot、facet 和未注册外部图当前不支持。

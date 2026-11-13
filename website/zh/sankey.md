## 先准备真实来源 {#input}

{{prepare:foundationdata}}

## 镜像双列关系 {#sankey}

Sankey 将每个来源在左右两列各显示一次。左侧 query 端点连接右侧 subject 端点，
互反输入行保持独立；默认排除同实体自身连接，`self_links=TRUE` 可启用。
宽度由两端真实区间决定，表达的是关系，不是流量守恒。本版由坐标系管理左右副本，
不接受 `between`，也不支持共享分支。

{{code:foundationsankey}}

{{figure:foundationsankey}}

## 只改变列间距 {#column-gap}

{{code:coursesankeygap}}

![仅改变左右列间距](../assets/coursesankeygap.png)

## 核对来源与实例 {#instances}

```r
get_chord_layout(p_foundationsankey)$display_instances
```

两条来源在左右各显示一次；不能手工复制 seq_data 冒充新观测。互反记录保持独立。self_links = TRUE 只启用已有自身记录，不补造关系。seq_order 每个来源只能出现一次；between 和 shared branches 不支持。

练习：从一条 ribbon 的左右端找回 query 和 subject，并解释为什么两端总宽不需守恒。

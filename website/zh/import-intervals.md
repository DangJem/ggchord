{{code:courseblockao}}

{{prepare:setup}}

先解压[真实输入包](../assets/tutorial-inputs.zip)，工作目录须含 examples/；每节单列可选依赖。

## BED 和 PAF：坐标只转换一次 {#zero-based}
BED 的起点从 0 开始，终点不包含在区间内。下例读入一条 pUC19 复制起点记录，比较读入前后的坐标。

{{code:formats}}
{{output:formats}}

PAF 提供两端全长和区间。先选出 mapq ≥ 30 的比对，再由两端长度生成骨架表。

{{code:formatpaf}}
{{output:formatpaf}}

{{code:drawpaf}}
{{figure:paf}}

反向比对用递减的 subject 端点表达。绘图使用读取器转换后的坐标。
## 可选的区间对象 {#ranges}

第一段把序列 a 上的 10–30 bp 转成特征表。第二段把 a 上的 10–30 bp 与 b 上的 40–60 bp 连成一条记录。GRanges 表示一个区间，GInteractions 表示一对区间。先安装 GenomicRanges 和 InteractionSet，完整序列长度由 seqlengths 提供。

{{code:frontimportgr}}

{{output:frontimportgr}}

{{code:frontimportgi}}

{{output:frontimportgi}}
## BED双端区间 {#bedpe}

`read_bedpe()`读取六个必填的制表符分隔字段，以及可选name、score、两条链和额外元数据。第一端保持query、第二端保持subject；score保留为字符串，不推断Identity或测序支持度。格式依据[BEDPE说明](https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format)。

下面将随包已观测比对的区间转写为BEDPE，用于坐标往返验证，不将它们当作配对测序证据。六列输入没有链证据，因此导入后的锚点链和相对方向保持未知。有真实链信息时才提供对应列，不由坐标升序猜测链。

{{code:courseblockap}}

负链锚点采用降序canonical端点，与配对区间适配器一致；同时保留原始start/end和链字段。name不视为唯一ID。额外列在`bedpe_extra`中保留，按`column_11`起命名。未知染色体（`.`）或边界（`-1`）默认报错；显式`unmapped="keep"`为审计保留不完整行，适配或绘图前须筛选`q_mapped & s_mapped`。非法已知坐标和零宽区间始终报错。区间Ribbon没有Identity时，显式设置固定填充，不制造pident。

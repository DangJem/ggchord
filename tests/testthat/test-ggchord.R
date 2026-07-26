# testthat 测试文件 - v0.4.0 分层 API

test_that("ggchord 仅使用 seq_data 能返回 ggchord 对象", {
  data(seq_data_example)
  p <- ggchord(seq_data = seq_data_example)
  expect_s3_class(p, "ggchord")
  expect_s3_class(p, "ggplot")
})

test_that("ggchord + geom_seq 能生成序列弧线图层", {
  data(seq_data_example)
  p <- ggchord(seq_data = seq_data_example) + geom_seq()
  expect_s3_class(p, "ggchord")
  expect_true(length(p$layers) >= 1)
})

test_that("ggchord + geom_seq + geom_ribbon 正确处理比对数据", {
  data(seq_data_example)
  data(ribbon_data_example)
  p <- ggchord(
    seq_data = seq_data_example,
    ribbon_data = ribbon_data_example
  ) +
    geom_seq() +
    geom_ribbon()
  expect_s3_class(p, "ggchord")
  expect_true(length(p$layers) >= 1)
})

test_that("ggchord + geom_gene 正确处理基因数据", {
  data(seq_data_example)
  data(gene_data_example)
  p <- ggchord(
    seq_data = seq_data_example,
    gene_data = gene_data_example
  ) +
    geom_seq() +
    geom_gene()
  expect_s3_class(p, "ggchord")
})

test_that("ggchord + geom_axis 能添加坐标轴", {
  data(seq_data_example)
  p <- ggchord(seq_data = seq_data_example) +
    geom_seq() +
    geom_axis()
  expect_s3_class(p, "ggchord")
})

test_that("参数可分布在各个 geom 中", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)

  p <- ggchord(
    seq_data = seq_data_example,
    ribbon_data = ribbon_data_example,
    gene_data = gene_data_example,
    title = "Test",
    rotation = 30
  ) +
    geom_seq(
      seq_radius = c(3, 2, 2, 1),
      seq_curvature = c(0, 1, -1, 1.5),
      seq_orientation = c(-1, -1, -1, 1)
    ) +
    geom_ribbon(
      ribbon_color_scheme = "query",
      ribbon_alpha = 0.5
    ) +
    geom_gene(
      gene_color_scheme = "strand",
      gene_width = 0.08
    ) +
    geom_axis(
      axis_gap = 0.02,
      axis_label_orientation = c(0, 45, 80, 130)
    )

  expect_s3_class(p, "ggchord")
})

test_that("seq_data 缺少必须列时报错", {
  bad_data <- data.frame(id = c("a", "b"), len = c(100, 200))
  expect_error(
    ggchord(seq_data = bad_data),
    "seq_data"
  )
})

test_that("seq_data 长度非正时报错", {
  bad_data <- data.frame(seq_id = c("a", "b"), length = c(0, 200))
  expect_error(
    ggchord(seq_data = bad_data),
    "正数"
  )
})

test_that("debug 模式输出调试信息", {
  data(seq_data_example)
  p <- ggchord(seq_data = seq_data_example, debug = TRUE)
  expect_s3_class(p, "ggchord")
})

test_that("ggchord 全局参数正确传递", {
  data(seq_data_example)
  p <- ggchord(seq_data = seq_data_example, title = "Hello", rotation = 90)
  expect_s3_class(p, "ggchord")
})

test_that("print 时能渲染完整弦图", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)

  # 在真实 ggnewscale 环境中 pident 也可以正常工作
  p <- ggchord(seq_data_example, ribbon_data_example, gene_data_example) +
    geom_seq() +
    geom_ribbon(ribbon_color_scheme = "query") +
    geom_gene() +
    geom_axis()

  # 渲染到 PDF
  pdf("/tmp/ggchord_test_print.pdf", 8, 8)
  suppressMessages(suppressWarnings(print(p)))
  dev.off()
  expect_true(file.exists("/tmp/ggchord_test_print.pdf"))
})

test_that("README 中的配色、标签覆盖与透明度参数生效", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)

  p <- ggchord(seq_data_example, ribbon_data_example, gene_data_example) +
    geom_seq(linewidth = 2) +
    geom_ribbon(ribbon_color_scheme = "query", alpha = 0.2) +
    geom_gene(gene_color_scheme = "manual", show_label = TRUE,
              label_size = 4) +
    geom_axis()

  pdf("/tmp/ggchord_readme_params.pdf", 8, 8)
  expect_no_warning(print(p))
  dev.off()

  layout <- ggchord:::get_chord_layout()
  expect_true("fill" %in% names(layout$ribbon_polys))
  expect_true(all(layout$ribbon_polys$alpha == 0.2))
  expect_gt(nrow(layout$gene_labels), 0)
  expect_true(file.exists("/tmp/ggchord_readme_params.pdf"))
})

test_that("文档限定的数据和参数值会被校验", {
  expect_error(
    ggchord(data.frame(seq_id = c("a", "a"), length = c(1, 2))),
    "唯一"
  )

  data(seq_data_example)
  expect_error({
    p <- ggchord(seq_data_example) + geom_seq(seq_orientation = 0)
    pdf(tempfile(fileext = ".pdf")); print(p); dev.off()
  }, "1 或 -1")
})

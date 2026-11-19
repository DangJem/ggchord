# Use the same topic memberships for the index and native Rd pages.
reference_sections <- yaml::read_yaml('_pkgdown.yml')$reference
reference_sections <- lapply(reference_sections, function(section) {
  section$topics <- gsub('`', '', section$contents, fixed = TRUE)
  section$files <- unname(alias_topics[section$topics])
  stopifnot(!anyNA(section$files))
  section
})
reference_files <- unlist(lapply(reference_sections, '[[', 'files'), use.names = FALSE)
stopifnot(!anyDuplicated(reference_files),
          setequal(reference_files, tools::file_path_sans_ext(basename(topic_files))))

# Parameter names are scoped to the function whose public result was checked.
flexible_parameters <- jsonlite::fromJSON("website/flexible-parameters.json",
                                          simplifyVector = FALSE)

mark_flexible_parameters <- function(main, topic, lang) {
  parameters <- flexible_parameters[[topic]]

  zh <- identical(lang, "zh")
  badge <- if (zh) "【支持灵活参数】" else "[Flexible parameter]"
  anchors <- c(sequence = "scalar", strand = "width", `strand-offset` = "strand",
               `fixed-label` = "scope", band = "scope")
  found <- character()
  for (term in xml2::xml_find_all(main, ".//dt[starts-with(@id, 'arg-')]")) {
    arguments <- trimws(strsplit(trimws(xml2::xml_text(term)), ",", fixed = TRUE)[[1]])
    permalink <- xml2::xml_find_first(term, "./a[@class='anchor']")
    xml2::xml_remove(xml2::xml_contents(term))
    for (i in seq_along(arguments)) {
      argument <- arguments[i]
      if (i > 1L) xml2::xml_add_child(term, "span", ", ")
      name <- xml2::xml_add_child(term, "span", class = "argument-term")
      xml2::xml_add_child(name, "code", argument, id = paste0("parameter-", argument),
                          class = "argument-anchor")
      kind <- parameters[[argument]]
      if (!is.null(kind)) {
        stopifnot(kind %in% names(anchors))
        xml2::xml_add_child(name, "a", badge, class = "flexible-parameter",
                            href = paste0("../parameters.html#", anchors[[kind]]))
        found <- c(found, argument)
      }
    }
    if (!inherits(permalink, "xml_missing")) xml2::xml_add_child(term, permalink)
  }
  if (!setequal(found, names(parameters))) {
    stop("Missing documented arguments in ", topic, ": ",
         paste(setdiff(names(parameters), found), collapse = ", "))
  }
}

add_reference_navigation <- function(doc, main, slug, lang) {
  zh <- lang == 'zh'
  label <- if (zh) 'API 参考' else 'API reference'
  if (slug == 'index') {
    xml2::xml_set_text(xml2::xml_find_first(main, './/h1'), label)
    xml2::xml_set_text(xml2::xml_find_first(doc, '//title'), paste(label, '• ggchord'))
    headings <- xml2::xml_find_all(main, './div/h2')
    stopifnot(length(headings) == length(reference_sections))
    for (i in seq_along(headings)) {
      section <- reference_sections[[i]]
      heading <- headings[[i]]
      xml2::xml_set_attr(heading, 'id', section$id)
      xml2::xml_set_attr(xml2::xml_find_first(heading, './a'), 'href', paste0('#', section$id))
      if (!is.null(section$previous_id)) {
        xml2::xml_add_sibling(heading, 'span', id = section$previous_id, .where = 'before')
      }
      wrapper <- xml2::xml_parent(heading)
      xml2::xml_set_attr(wrapper, 'class', 'section level2 reference-group')
      # pkgdown emits a title row and a subsequent topic row for each section.
      topics <- xml2::xml_find_first(wrapper, 'following-sibling::div[1]/dl')
      stopifnot(length(xml2::xml_find_all(topics, './dt')) == length(section$files))
      xml2::xml_remove(xml2::xml_find_all(wrapper, './dl'))
      empty <- xml2::xml_parent(topics)
      xml2::xml_add_child(wrapper, topics)
      xml2::xml_remove(empty)
    }
    filter <- xml2::xml_add_sibling(xml2::xml_parent(xml2::xml_find_first(main, './/h1')),
                                     'div', class = 'reference-filter', hidden = 'hidden')
    field <- xml2::xml_add_child(filter, 'label', if (zh) '按函数名或用途筛选' else 'Filter by name or purpose')
    xml2::xml_set_attr(field, 'for', 'reference-filter')
    xml2::xml_add_child(filter, 'input', id = 'reference-filter', type = 'search', autocomplete = 'off')
    xml2::xml_add_child(filter, 'output', id = 'reference-filter-count', 'aria-live' = 'polite')
    return(invisible())
  }
  mark_flexible_parameters(main, slug, lang)
  section <- Filter(function(section) slug %in% section$files, reference_sections)[[1L]]
  idx <- match(slug, section$files)
  nearby <- head(order(abs(seq_along(section$files) - idx))[-1L], 3L)
  more <- xml2::xml_add_child(main, 'section', class = 'reference-related')
  xml2::xml_add_child(more, 'h2', if (zh) '相关条目' else 'Related entries', id = 'ref-related')
  links <- xml2::xml_add_child(more, 'ul')
  for (i in sort(nearby)) {
    item <- xml2::xml_add_child(links, 'li')
    link <- xml2::xml_add_child(item, 'a', href = paste0(section$files[i], '.html'))
    xml2::xml_add_child(link, 'code', section$topics[i])
  }
  item <- xml2::xml_add_child(links, 'li')
  xml2::xml_add_child(item, 'a', if (zh) section$title_zh else section$title,
                      href = paste0('index.html#', section$id))
  if (slug == 'ggchord_api_contract') {
    item <- xml2::xml_add_child(links, 'li')
    xml2::xml_add_child(item, 'a', if (zh) '版本选择' else 'Version choice', href = '../migration.html')
  }
}

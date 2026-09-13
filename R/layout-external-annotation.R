# Perimeter annotation coordination for circular maps.

# Restriction-site labels already own a polar contour/fan solver that keeps
# sparse labels at their natural angles and opens dense clusters locally. Do
# not flatten those labels onto Cartesian side rails when feature callouts are
# also present: that destroys the genomic contour and produces very long
# leaders. Feature labels exhaust their internal feature and label tracks
# before becoming external callouts, and each layer-local solver retains
# ownership of its final geometry.
ggchord_share_external_annotations <- function(registry, layout) {
  registry
}

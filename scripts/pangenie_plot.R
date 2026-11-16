require(RIdeogram)
args = (commandArgs(TRUE))

add_axis_labels <- function(svg_file = "chromosome.svg") {
    svg <- readLines(svg_file, warn = FALSE)
    labels <- paste0(
        '<text x="50%" y="64%" text-anchor="middle" ',
        'font-family="Arial, sans-serif" font-size="32" font-weight="700" fill="#667085">',
        'Chromosome number</text>\n',
        '<text x="-33%" y="32" text-anchor="middle" transform="rotate(-90)" ',
        'font-family="Arial, sans-serif" font-size="32" font-weight="700" fill="#667085">',
        'Length</text>'
    )
    svg <- sub('</svg>', paste0(labels, '\n</svg>'), svg, fixed = TRUE)
    writeLines(svg, svg_file)
}

data(human_karyotype, package="RIdeogram")

contigs <- read.table("contigs.tsv", sep = "\t", header = T, stringsAsFactors = F)

SV <- read.table(args[1], sep = "\t", header = T, stringsAsFactors = F)

ideogram(karyotype = human_karyotype, overlaid = contigs, label = SV, label_type = "marker", colorset1 = c("#627A9D", "#77AAE4"))
add_axis_labels("chromosome.svg")

convertSVG("chromosome.svg", device = "png")

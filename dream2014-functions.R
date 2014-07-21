"read.gct.file" = function(filename) {
  df = read.delim(filename, skip = 2, stringsAsFactors = F)
  df
}

"read.annotations" = function(filename) {
  df = read.delim(filename, stringsAsFactors = F)
  df
  # for expression: columns are Probe_Name, Gene_Name
  # for copy number: columns are Gene_name, Num_chr, txStart, txEnd
}

"read.features.list" = function(filename) {
  l = readLines(filename)
  l
}

"write.res.gct" = function(out.mat, gene.names, cell.line.names, outfile = "output.txt") {
  sink(outfile)
  cat("#1.2\n")
  cat(ncol(out.mat), "\t",  nrow(out.mat), "\n")
  sink()
  mat.res = t(out.mat)
  colnames(mat.res) = cell.line.names
  df = data.frame(Name = gene.names, Description = gene.names, mat.res)
  write.table(df, file = outfile, append = T, sep = "\t", row.names = F, col.names = T, quote = F)
}


# constructor function to create a zi object
zi <- function(ziInput, ziModel, ziOutput, weights) {
  structure(list("ziInput" = ziInput, "ziModel" = ziModel, "ziOutput" = ziOutput, "weights" = weights), class = c("zi","matrix"))
}

ps_replaced <- function(result_zi)
{
  ps_new <- result_zi$ziInput
  new_otu <- otu_table(result_zi$ziOutput, taxa_are_rows(result_zi$ziInput))
  otu_table(ps_new)<-new_otu
  return(ps_new)
}


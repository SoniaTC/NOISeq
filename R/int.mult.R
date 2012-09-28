int.mult <- function(lista, todos = NULL) {
  
  if(is.null(todos)) {
    todos <- unlist(lista)
  }

  comunes <- todos

  for(i in 1:length(lista)) {
    comunes <- intersect(comunes, lista[[i]])
  }

  comunes
}

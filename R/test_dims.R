find_dimensions = function(number) {
  val1 <- floor(sqrt(number))
  val2 <- ceiling(sqrt(number))
  #res1 <- val1
  #res2 <- val2
  
  while (val1*val2 < number) {
    val2 <- val2 + 1
    if (val1 * val2 < number) {
      val1 <- val1 - 1
    }
  }
  
  return(c(val1, val2))
  
}
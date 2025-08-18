library(stringr)

longest_common_substring <- function(strings, min_length = 3) {
  if (length(strings) == 0) return("")
  
  # Start with the shortest string
  shortest <- strings[which.min(nchar(strings))]
  len <- nchar(shortest)
  
  # Binary search on substring length
  is_common <- function(sub_len) {
    substrings <- unique(sapply(1:(len - sub_len + 1), function(start) {
      substr(shortest, start, start + sub_len - 1)
    }))
    any(sapply(substrings, function(sub) all(grepl(sub, strings))))
  }
  
  # Perform binary search for the longest common substring
  left <- 0
  right <- len
  best_substring <- ""
  
  while (left <= right) {
    mid <- (left + right) %/% 2
    substrings <- unique(sapply(1:(len - mid + 1), function(start) {
      substr(shortest, start, start + mid - 1)
    }))
    
    valid_substrings <- substrings[sapply(substrings, function(sub) all(grepl(sub, strings)))]
    
    if (length(valid_substrings) > 0) {
      best_substring <- valid_substrings[which.max(nchar(valid_substrings))]  # Store the longest valid substring
      left <- mid + 1  # Try finding a longer one
    } else {
      right <- mid - 1  # Reduce the search space
    }
  }

  # Check if the found substring meets the minimum length requirement
  if (nchar(best_substring) < min_length) {
    best_substring <- ""
  }
  
  return(best_substring)
}
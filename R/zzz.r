# .onLoad <- function(libname, pkgname) {
#   memo.bic <<- memoise::memoise(memo.bic,
#                                 cache = cachem::cache_mem(evict = "lru"),
#                                 omit_args = c("param_start"))
#   IBIS.Z <<- memoise::memoise(IBIS.Z,cache = cachem::cache_mem(evict = "lru"),
#                               omit_args = c("sampl","current_set"))
# }

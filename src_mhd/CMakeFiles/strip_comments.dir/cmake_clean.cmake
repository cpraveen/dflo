FILE(REMOVE_RECURSE
  "CMakeFiles/strip_comments"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/strip_comments.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)

function(check_mindsplatter_replay_revision committed generated)
  set(_revision_pattern "msp-heavy-search-v[0-9]+")
  string(REGEX MATCH "${_revision_pattern}" _committed_revision "${committed}")
  string(REGEX MATCH "${_revision_pattern}" _generated_revision "${generated}")
  if(_committed_revision STREQUAL "" OR _generated_revision STREQUAL "")
    message(FATAL_ERROR
      "MindSplatter replay revision missing: committed ${_committed_revision}, "
      "generated ${_generated_revision}")
  endif()
  if(NOT _committed_revision STREQUAL _generated_revision)
    message(FATAL_ERROR
      "MindSplatter replay revision drift: committed ${_committed_revision}, "
      "generated ${_generated_revision}")
  endif()
  set(MINDSPLATTER_REPLAY_REVISION "${_generated_revision}" PARENT_SCOPE)
endfunction()

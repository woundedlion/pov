# Require the device inline-storage audit in core/animation/animation.h to name
# every concrete animation type. LARGEST_CONCRETE_ANIM_SIZE folds sizeof over a
# hand-written pack, and the per-type static_assert in Timeline::add_get only
# fires for types a given build actually add()s — so a new animation that no
# build in the tree instantiates is measured by nothing, which is the gap the
# audit exists to close.
#
# Animations are CRTP: `class X : public Base<X>`, rooted at AnimationBase.
# Each declaration is classified by the base's template argument:
#   * argument == the class name          -> concrete animation; must be packed.
#   * argument == the class name + '<...>' -> templated animation; its size
#                                             depends on the instantiation, so
#                                             the pack cannot name it.
#   * anything else (e.g. `Derived`)       -> intermediate CRTP base; its own
#                                             subclasses are animations too.
# The root set therefore grows to a fixpoint before the leaves are collected.
# -D args: CORE_DIR (path to core/), ANIM_HEADER (path to animation.h).

cmake_minimum_required(VERSION 3.29)

set(_roots "AnimationBase")
set(_concrete "")
set(_templated "")

file(GLOB_RECURSE _headers "${CORE_DIR}/*.h")
set(_decls "")
foreach(_hdr IN LISTS _headers)
  file(READ "${_hdr}" _text)
  # Strings first, so a `//` or `/*` inside a message cannot open a comment span.
  string(REGEX REPLACE "\"([^\"\\\\\n]|\\\\.)*\"" "\"\"" _text "${_text}")
  string(REGEX REPLACE "/\\*[^*]*\\*+([^/*][^*]*\\*+)*/" "\n" _text "${_text}")
  string(REGEX REPLACE "//[^\n]*" "" _text "${_text}")
  # clang-format wraps a long base-clause onto its own line, so allow newlines
  # around the colon and after `public`.
  string(REGEX MATCHALL
    "class[ \t\r\n]+[A-Za-z0-9_]+[ \t\r\n]*:[ \t\r\n]*public[ \t\r\n]+[A-Za-z0-9_]+[ \t\r\n]*<[^;{]*>"
    _matches "${_text}")
  list(APPEND _decls ${_matches})
endforeach()

if(NOT _decls)
  message(FATAL_ERROR
    "no CRTP class declarations found under ${CORE_DIR}; this check cannot "
    "verify the animation pack and must not report success.")
endif()

# Fixpoint: an intermediate base discovered on one pass admits its own
# subclasses on the next. Four passes cover any hierarchy depth in the tree.
foreach(_pass RANGE 3)
  foreach(_decl IN LISTS _decls)
    string(REGEX REPLACE "^class[ \t\r\n]+([A-Za-z0-9_]+).*" "\\1" _name "${_decl}")
    string(REGEX REPLACE "^class[ \t\r\n]+[A-Za-z0-9_]+[ \t\r\n]*:[ \t\r\n]*public[ \t\r\n]+([A-Za-z0-9_]+).*"
      "\\1" _base "${_decl}")
    if(NOT _base IN_LIST _roots)
      continue()
    endif()
    # Base template argument: everything between the outermost <> of the base.
    string(REGEX REPLACE "^[^<]*<(.*)>$" "\\1" _arg "${_decl}")
    string(STRIP "${_arg}" _arg)
    if(_arg STREQUAL "${_name}")
      if(NOT _name IN_LIST _concrete)
        list(APPEND _concrete "${_name}")
      endif()
    elseif(_arg MATCHES "^${_name}[ \t\r\n]*<")
      if(NOT _name IN_LIST _templated)
        list(APPEND _templated "${_name}")
      endif()
    elseif(NOT _name IN_LIST _roots)
      list(APPEND _roots "${_name}")
    endif()
  endforeach()
endforeach()

# The pack: the argument list of largest_sizeof<...>() feeding
# LARGEST_CONCRETE_ANIM_SIZE.
file(READ "${ANIM_HEADER}" _anim)
string(REGEX REPLACE "/\\*[^*]*\\*+([^/*][^*]*\\*+)*/" "\n" _anim "${_anim}")
string(REGEX REPLACE "//[^\n]*" "" _anim "${_anim}")
if(NOT _anim MATCHES "LARGEST_CONCRETE_ANIM_SIZE[ \t\r\n]*=[ \t\r\n]*largest_sizeof<([^>]*)>")
  message(FATAL_ERROR
    "cannot find `LARGEST_CONCRETE_ANIM_SIZE = largest_sizeof<...>` in "
    "${ANIM_HEADER}; the inline-storage audit this check pins is gone or was "
    "renamed.")
endif()
string(REPLACE "Animation::" "" _pack_text "${CMAKE_MATCH_1}")
string(REGEX MATCHALL "[A-Za-z0-9_]+" _packed "${_pack_text}")

set(_missing "")
foreach(_name IN LISTS _concrete)
  if(NOT _name IN_LIST _packed)
    list(APPEND _missing "${_name}")
  endif()
endforeach()

set(_stale "")
foreach(_name IN LISTS _packed)
  if(NOT _name IN_LIST _concrete)
    list(APPEND _stale "${_name}")
  endif()
endforeach()

if(_missing)
  string(REPLACE ";" ", " _missing_list "${_missing}")
  message(FATAL_ERROR
    "concrete animation type missing from LARGEST_CONCRETE_ANIM_SIZE: "
    "${_missing_list}. Every non-templated AnimationBase subclass must be named "
    "in the pack in ${ANIM_HEADER}, or its size is never checked against the "
    "TimelineEvent inline-storage budget in a build that does not add() it.")
endif()

if(_stale)
  string(REPLACE ";" ", " _stale_list "${_stale}")
  message(FATAL_ERROR
    "LARGEST_CONCRETE_ANIM_SIZE names a type that is not a concrete "
    "AnimationBase subclass: ${_stale_list}. Remove it from the pack in "
    "${ANIM_HEADER}, or the audit measures a type the timeline never stores.")
endif()

list(LENGTH _concrete _nconcrete)
list(LENGTH _templated _ntemplated)
message(STATUS
  "animation pack: all ${_nconcrete} concrete AnimationBase subclasses are "
  "named in LARGEST_CONCRETE_ANIM_SIZE (${_ntemplated} templated ones are "
  "sized per instantiation and excluded)")

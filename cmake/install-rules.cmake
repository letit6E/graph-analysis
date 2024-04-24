install(
    TARGETS graph-analysis_exe
    RUNTIME COMPONENT graph-analysis_Runtime
)

if(PROJECT_IS_TOP_LEVEL)
  include(CPack)
endif()

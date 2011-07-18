###############
# main target #
###############

set (sources ${src_data} ${src_linesearch} ${src_loss} ${src_model} ${src_regularizer} ${src_solver} ${src_tools} ${src_utilities})

#add_executable(linear-bmrm-predict linear-bmrm-predict.cpp ${sources})
#target_link_libraries(linear-bmrm-predict ${CORELIBS})

add_executable(linear-bmrm-train linear-bmrm-train.cpp ${sources})
target_link_libraries(linear-bmrm-train ${CORELIBS})

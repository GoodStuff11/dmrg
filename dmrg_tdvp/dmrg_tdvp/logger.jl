function set_logger(output_path)
    logger = FormatLogger(open(joinpath(output_path,"log"), "a+")) do io, args
        print(io, "Time: ")
        Dates.format(io, now(), dateformat"yyyy-mm-dd HH:MM:SS")
        print(io, " | ")
        print(io, args.level, ": ")
        println(io, args.message)
    end

    global_logger(logger)
end


# mutable struct Logger
#     path::String
#     filename::String

#     function Logger(path::String)
#         new(path, "log")
#     end
#     function Logger(path::String, filename::String)
#         new(path, filename)
#     end
# end

# function Base.:<(log::Logger, strings::Vector{Any})
#     file = open(joinpath(log.path, log.filename),"a")
#     for s in strings
#         if s isa Tuple
#             println(file, s...)
#         else
#             println(file,s)
#         end
#     end
#     close(file)
# end

# function Base.:<(log::Logger, string::Union{Any, Tuple{Any}})
#     file = open(joinpath(log.path, log.filename),"a")
#     if string isa Tuple
#         println(file,string...)
#     else
#         println(file,string)
#     end
#     close(file)
# end

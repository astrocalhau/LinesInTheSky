function write_fits(filename, results)
    data = OrderedDict{String, Vector}()
    for cname in names(results)
        col = results[:, cname]
        typ = nonmissingtype(eltype(col))
        if typ <: AbstractFloat
            data[cname] = replace(col, missing => NaN)
        elseif typ <: Integer
            data[cname] = replace(col, missing => -1)
        elseif typ <: String
            data[cname] = replace(col, missing => "")
        else
            error("Unsupported data type: $(typ)")
        end
    end
    f = FITS(filename, "w")
    write(f, data)
    close(f)
end

function run_analysis(input_path, output_path, catalog)
    # bash -c 'echo "`ls results/JSON | wc -w` / `ls input/ | wc -w`" | bc -l'
    SENTINEL = -1
    channel = Channel{Int64}(30)
    function consumer(ith::Int, channel::Channel)
        while true
            i = take!(channel)
            if i == SENTINEL
                put!(channel, SENTINEL) # tell other threads to quit
                break                   # quit this thread
            end
            try
                analyze_spec(input_path, output_path, catalog[i, :])
            catch err
                display(err)
                println("Failed to fit spectrum for i=$i")
            end
        end
    end
    tasks = [@spawn consumer(ith, channel) for ith in 1:nthreads()]
    put!.(Ref(channel), 1:nrow(catalog))  # tell threads which row to analyze
    put!(     channel , SENTINEL)         # tell threads to quit
    wait.(tasks)                          # wait for threads to terminate
end

function read_all_results(output_path, catalog)
    SENTINEL = -1
    channel = Channel{Int64}(30)

    out = [DataFrame() for i in 1:(Threads.nthreads(:interactive) .+ Threads.nthreads(:default))]
    function consumer(ith::Int, channel::Channel)
        df = out[ith]
        while true
            i = take!(channel)
            if i == SENTINEL
                put!(channel, SENTINEL) # tell other threads to quit
                break                   # quit this thread
            end
            try
                cur = read_results(output_path, catalog[i, :])
                append!(df, cur, cols=:union)
            catch err
                display(err)
                println("Error reading results for i=$i")
            end
        end
    end

    tasks = [@spawn consumer(ith, channel) for ith in 1:nthreads()]
    put!.(Ref(channel), 1:nrow(catalog))  # tell threads which row to analyze
    put!(     channel , SENTINEL)         # tell threads to quit
    wait.(tasks)                          # wait for threads to terminate
    return vcat(out..., cols=:union)
end

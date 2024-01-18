@testset "__init__ env vars" begin
    showtestset()

    upstream_env_vars = [
        "DATAPATH",
        "CYCLEFOLD_DATAPATH",
    ]
    our_env_vars = ["RNASTRUCTURE_JL_" * e for e in upstream_env_vars]
    env_var_mapping = Dict(e => "RNASTRUCTURE_JL_" * e for e in upstream_env_vars)
    all_env_vars = [upstream_env_vars..., our_env_vars...]

    function delete_all_env_vars()
        for e in all_env_vars
            delete!(ENV, e);
        end
    end

    # save relevant ENV vars so that we can change them here for
    # testing and then restore them later
    saved_env_vars = Dict{String,String}()
    for e in all_env_vars
        if haskey(ENV, e)
            saved_env_vars[e] = ENV[e]
            delete!(ENV, e)
        end
    end

    @test RNAstructure.__init__() == nothing
    delete_all_env_vars()

    # setting DATAPATH, etc
    for e in upstream_env_vars
        ENV[e] = "non-existent-dir-path"
        warn_msg = ("RNAstructure: $e env var set, replacing with $(saved_env_vars[e])\n"
                    * "To override $e used by RNAstructure, set the RNASTRUCTURE_JL_$e env var")
        @test (@test_logs (:warn, (warn_msg)) RNAstructure.__init__()) == nothing
        delete_all_env_vars()
    end

    # setting RNASTRUCTURE_JL_*
    for (upstream_env, our_env) in env_var_mapping
        ENV[our_env] = "non-existent-dir-path"
        info_msg = ("Setting ENV[\"$upstream_env\"] = ENV[\"$our_env\"]")
        @test (@test_logs (:info, (info_msg)) RNAstructure.__init__()) == nothing
        delete_all_env_vars()
    end

    # setting both upstream_env_vars and our_env_vars, our_env_vars
    # should have precedence
    for (upstream_env, our_env) in env_var_mapping
        ENV[upstream_env] = "non-existent-dir-path"
        ENV[our_env] = "non-existent-dir-path"
        info_msg = ("Setting ENV[\"$upstream_env\"] = ENV[\"$our_env\"]")
        @test (@test_logs (:info, (info_msg)) RNAstructure.__init__()) == nothing
        delete_all_env_vars()
    end

    # restore env vars that we previously unset
    for e in all_env_vars
        if haskey(saved_env_vars, e)
            ENV[e] = saved_env_vars[e]
        end
    end
end

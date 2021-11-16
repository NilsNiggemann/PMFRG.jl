function capture_stdout(f)
    #https://discourse.julialang.org/t/consistent-way-to-suppress-solver-output/20437/6
    stdout_orig = stdout
    (rd, wr) = redirect_stdout()
    f()
    close(wr)
    redirect_stdout(stdout_orig)
    read(rd, String)
end

function __precompile__()
    System = getPolymer(4)
    Par = Params(System = System,N=16,T=1.2,MinimalOutput=false,usesymmetry = true,accuracy= 1E-4,Lam_min = 10)
    mainF = UniqueDirName("temp")*".h5"
    CheckF = UniqueDirName("temp_ch")
    Solution,saved_values = SolveFRG(Par,OneLoop(),MainFile = mainF,CheckpointDirectory = CheckF,method = DP5(),VertexCheckpoints = [],overwrite_Checkpoints=true)
    rm(mainF)
    Solution,saved_values = SolveFRG(Par,TwoLoop(),MainFile = mainF,CheckpointDirectory = CheckF,method = DP5(),VertexCheckpoints = [],overwrite_Checkpoints=true)
    rm.((mainF,CheckF),recursive=true)
    # Solution,saved_values = SolveFRG_Checkpoint("/storage/niggeni/test/Cubic_NLen=4_N=24_T=1.200000000000000/CurrentState.h5",getCubic,OneLoop(),MainFile="/storage/niggeni/test/main.h5",method = BS3(),VertexCheckpoints = [],overwrite_Checkpoints = true)
    return 
end

function __precompile__quiet__()
    capture_stdout(__precompile__)
    return
end
function networksetup(c;digits=3,nSC=2,nFC=1,N=2,normalise=false)
  
        SC,dist= get_paul_data_all(;nSC=1,nFC=1,type="control",ROI=N)
       
        N = size(SC,1)
        lags = dist./c
        return SC,dist,lags,N
   
    
end


function get_paul_data_all(;nSC=1,nFC=1,type="control",ROI=140)

    #FUNCTDIR="$DATADIR/Functional/$ROIDIR"
    
    
    SC = get_stuct_data(;n=nSC,ROI=ROI)
    dist = get_dist_data(;n=nSC,ROI=ROI)

    
    return SC,dist
end
    
    
function get_stuct_data(;n=1,ROI=140)

    if ROI in [18,64,140,246,503,673] 
        HOMEDIR = homedir()
        ROIDIR = "ROI"*"$ROI"
        DATADIR = "$HOMEDIR/PaulJuliaData_ALL"
        STRUCTDIR="$DATADIR/Structural/$ROIDIR"
        StructMats = readdir(STRUCTDIR)
        numMats = size(StructMats,1)
        if n > numMats
            @error "There are only $numMats structural matrices in this directory. Please choose n < 5"
        else

            SC = load("$STRUCTDIR/$(StructMats[n])","$(split(StructMats[n],".")[1])")

            SC = log.(SC)
            SC[SC .== -Inf] .= 0 
            SC = SC/maximum(SC)

            

            return SC
        end
    else
        @error "no such ROI size, choose an ROI in [18,64,140,246,503,673]"
    end

end

function get_dist_data(;n=1,ROI=140)
    if ROI in [18,64,140,246,503,673] 
        HOMEDIR = homedir()
        ROIDIR = "ROI"*"$ROI"
        DATADIR = "$HOMEDIR/PaulJuliaData_ALL"
        DISTDIR="$DATADIR/Distance/$ROIDIR"
        DistMats = readdir(DISTDIR)
        numMats = size(DistMats,1)
        if n > numMats
            @error "There are only $numMats distance matrices in this directory. Please choose n < 5"
        else

            dist = load("$DISTDIR/$(DistMats[n])","$(split(DistMats[n],".")[1])")

            return dist
        end
    else
        @error "no such ROI size, choose an ROI in [18,64,140,246,503,673]"
    end
end

function get_functonal_data(;n=1,type="control",ROI=140)
    if ROI in [18,64,140,246,503,673] 
        HOMEDIR = homedir()
        ROIDIR = "ROI"*"$ROI"
        DATADIR = "$HOMEDIR/PaulJuliaData_ALL/Functional"
        if lowercase(type) == "control"
            FUNCTIONALDIR="$DATADIR/$ROIDIR/Control_Groups"
        elseif lowercase(type) == "stimulated"
            FUNCTIONALDIR="$DATADIR/$ROIDIR/Stimulated_Groups"
        end

        FunctionalDIR = readdir(FUNCTIONALDIR)[n]
        
        
        #println(FUNCTIONALDIR*"/"*FunctionalDIR)
        
        FunctionalMats = readdir(FUNCTIONALDIR*"/"*FunctionalDIR)
        FunctionalMats = FunctionalMats[FunctionalMats .!== "MISSING_ROIs.jld"]
        missingROIs = load(FUNCTIONALDIR*"/"*FunctionalDIR*"/"*"MISSING_ROIs.jld","MISSING_ROIs")

        sizeFC = ROI - size(missingROIs,1)
        FC = zeros(sizeFC,sizeFC)

        numMats = size(FunctionalMats,1)
  
        for i in FunctionalMats
            FC += load("$FUNCTIONALDIR/$FunctionalDIR/$i","$(split(i,".")[1])")
        
        end


        return FC./numMats,missingROIs
    else
        @error "no such ROI size, choose an ROI in [18,64,140,246,503,673]"
    end
end
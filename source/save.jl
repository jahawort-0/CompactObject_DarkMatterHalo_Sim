module Save
    using CSV, Dates
    
    """Formats the filename as YYYYMMDD_M_WD=*M_WD*_M_NS=*M_NS*_A=*A*.csv
		if a filename is not provided"""
	function format_filename(M_WD::Float64,M_NS::Float64, A::Float64, mode; filename::Union{String,Int}=-1)::String
		#
        M_WD=round(M_WD,digits=3)
        M_NS=round(M_NS,digits=3)
        A=round(A,digits=3)
        if mode==:isotropic
            m_str="Iso"
        elseif mode==:jeans
             m_str="Jeans"
        elseif mode==:circumbinary
            m_str="CR"
        end
            
		if filename==-1 
		        today_str = Dates.format(Dates.now(), "yyyymmdd")
		        base_filename = "Outputs_3/$(today_str)_M_WD=$(M_WD)_M_NS=$(M_NS)_A=$(A)_mode=$(m_str).csv"
		        filename = base_filename
		        suffix = 2
		        while isfile(filename)
		            filename = "Outputs_3/$(today_str)_M_WD=$(M_WD)_M_NS=$(M_NS)_A=$(A)_$(suffix).csv"
		            suffix += 1
				end
		end
		return filename
	end
		
	"""Writes array of format given in Integrate.integrate() to a file."""
	function save_as_csv(datumses::AbstractMatrix{<:Real},
					  M_WD::Float64,M_NS::Float64,A::Float64,mode; filename::Union{String,Int}=-1)::Nothing
		
		
		filename=format_filename(M_WD,M_NS,A,mode;filename=filename)
		# Writes the file to a CSV titled filename

		isdir("Outputs_3") || mkdir("Outputs_3")
		#creates output directory if not already extant
		
			CSV.write(filename, 
						(t = datumses[:,1],
     					ddI_xx = datumses[:,2],
     					ddI_xy = datumses[:,3],
						a = datumses[:,4],
     					M_WD = datumses[:,5],
     					M_NS = datumses[:,6],
						R_WD = datumses[:,7],
     					J = datumses[:,8],
     					θ = datumses[:,9]))
		return nothing
	end

    function save_as_csv_append(datumses::AbstractMatrix{<:Real},
                            M_WD::Float64, M_NS::Float64, A::Float64, mode;
                            filename::Union{String,Int} = -1)::Nothing

        # Construct filename only if needed
        if filename isa Int || filename == ""
            filename = Save.format_filename(M_WD, M_NS, A, mode; filename=filename)
        end
    
        isdir("Outputs_3") || mkdir("Outputs_3")
    
        header_needed = !isfile(filename)
    
        open(filename, "a") do io
            # Write header manually if needed
            if header_needed
                println(io, "t,ddI_xx,ddI_xy,a,M_WD,M_NS,R_WD,J,theta")
            end
    
            # Append each row individually
            for k in 1:size(datumses, 1)
                println(io, join(datumses[k, :], ","))
            end
        end
    
        return nothing
    end
end
function overlap_bp(as, ae)
	--print("Inside overlap_bp")
	--print(start)
	--print(stop)
	--print(as)
	--print(ae)
	local start = start+1
	local stop = stop -1
	local diff = math.min(ae,stop) - math.max(as,start)
	diff = diff+1
	--print(result)
	--print("End of overlap_bp\n")
	return string.format("%.4f",diff)
end

function query_overlap(qas, qae)
	--print("Inside query_overlap")
	local qstart = start+1
	local qstop = stop -1
	--
       	diff = math.min(qae,qstop) - math.max(qas,qstart)+1
       	query_ov = diff/((qstop-qstart)+1)
	--print(result)
	--print(query_ov)
	--print("\n")
	return string.format("%.4f",query_ov)
end

function db_overlap(das, dae)
	--print("Inside db_overlap")
	local qstart = start+1
	local qstop = stop-1
       	diff = math.min(dae,qstop) - math.max(das,qstart)+1
	db_ov = diff/((dae-das)+1)
	--print(result)
	--print(db_ov)
	--print("\n")
	return string.format("%.4f",db_ov)
end

function printOps(ss,se,flag)
	--print("Inside printOps")
	--print(ss)
	--print(se)
	local ss_t = {}
	local se_t = {}

	index = 1
	for ele in string.gmatch(ss,"[^,]+") do
		ss_t[index] = ele
		index = index+1
	end

	index = 1
	for ele in string.gmatch(se,"[^,]+") do
		se_t[index] = ele
		index = index+1
	end
	
	local result = {}
	local ovlap_query = {}
	local ovlap_db = {}

	for i =1,table.getn(ss_t) do 
		local sst = tonumber(ss_t[i])
		local set = tonumber(se_t[i])
		
		if flag==0 then
			ovlap = overlap_bp(sst, set)
		elseif flag==1 then
			ovlap = query_overlap(sst, set)
		elseif flag==2 then
			ovlap = db_overlap(sst, set)
		end
		result[i] = ovlap
		--overlap_bp[i] = ovlap
	end
	
	return table.concat(result,",")

end

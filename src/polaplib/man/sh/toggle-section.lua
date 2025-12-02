function Pandoc(doc)
	local showDetails = true

	-- Read metadata key
	if doc.meta.showDetails ~= nil then
		showDetails = doc.meta.showDetails == true
		io.stderr:write("showDetails set to: ", tostring(showDetails), "\n")
	end

	local new_blocks = {}

	for _, block in ipairs(doc.blocks) do
		if block.t == "Div" and block.classes:includes("toggle-section") then
			io.stderr:write("Found a div with class 'toggle-section'\n")
			if showDetails then
				io.stderr:write("Keeping this div\n")
				table.insert(new_blocks, block)
			else
				io.stderr:write("Removing this div\n")
				-- skip
			end
		else
			table.insert(new_blocks, block)
		end
	end

	return pandoc.Pandoc(new_blocks, doc.meta)
end

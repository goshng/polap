function include_file(el)
	-- Convert paragraph to plain text
	local text = pandoc.utils.stringify(el)

	-- Check if the paragraph matches "!include filename.md"
	local filename = text:match("^!include%s+(.+)$")

	if filename then
		local file = io.open(filename, "r")
		if file then
			local content = file:read("*all")
			file:close()

			-- Convert the included content into Pandoc document blocks
			return pandoc.read(content, "markdown").blocks
		else
			return {} -- Return an empty block if the file is missing
		end
	end
end

return {
	{ Para = include_file },
}

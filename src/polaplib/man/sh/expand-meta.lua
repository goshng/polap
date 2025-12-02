local meta = {}

function Meta(m)
	meta = m
	return m
end

function walk_block(block)
	-- Convert the block to plain text
	local raw = pandoc.utils.stringify(block)
	-- Replace all {{key}} with metadata
	local replaced = raw:gsub("{{(.-)}}", function(key)
		return meta[key] and pandoc.utils.stringify(meta[key]) or "{{" .. key .. "}}"
	end)
	-- Re-parse replaced content into block
	return pandoc.read(replaced).blocks[1]
end

return {
	{ Block = walk_block },
}

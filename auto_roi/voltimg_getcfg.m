function val = voltimg_getcfg(cfg, field, default)
%VOLTIMG_GETCFG Read a config field with a default fallback.
%   val = VOLTIMG_GETCFG(cfg, field, default) returns cfg.(field) if cfg is a
%   struct containing that field (and it is nonempty), otherwise default.

if isstruct(cfg) && isfield(cfg, field) && ~isempty(cfg.(field))
    val = cfg.(field);
else
    val = default;
end
end

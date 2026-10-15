function h=p44_sha256(value,isFile)
if nargin<2; isFile=false; end
md=java.security.MessageDigest.getInstance('SHA-256');
if isFile
    fid=fopen(value,'rb'); assert(fid>0,'Cannot read %s.',value);
    cleanup=onCleanup(@()fclose(fid));
    while true
        bytes=fread(fid,1024*1024,'*uint8'); if isempty(bytes); break; end
        md.update(bytes);
    end
else
    md.update(unicode2native(char(value),'UTF-8'));
end
raw=typecast(md.digest(),'uint8'); h=lower(reshape(dec2hex(raw,2).',1,[]));
end

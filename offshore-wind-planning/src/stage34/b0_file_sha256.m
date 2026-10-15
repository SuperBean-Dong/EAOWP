function digest = b0_file_sha256(path)
% Hash a local file with or without the JVM; paths are never evaluated as code.
assert(isfile(path),'Missing file: %s',path);
if usejava('jvm')
    md=java.security.MessageDigest.getInstance('SHA-256');
    fid=fopen(path,'r'); assert(fid>0); cleanup=onCleanup(@()fclose(fid));
    while ~feof(fid)
        bytes=fread(fid,4*1024*1024,'*uint8');
        if ~isempty(bytes); md.update(typecast(bytes,'int8')); end
    end
    value=typecast(md.digest(),'uint8');
    digest=lower(reshape(dec2hex(value,2).',1,[]));
else
    assert(ispc,'Without JVM the packaged hash helper requires Windows.');
    assert(~contains(path,'"')&&~contains(path,'%')&&~contains(path,newline));
    literalPath=strrep(char(path),'''','''''');
    command=sprintf(['powershell.exe -NoProfile -NonInteractive -Command "', ...
        '$ErrorActionPreference=''Stop''; $env:PSModulePath=$PSHOME+''/Modules''; ', ...
        '(Get-FileHash -LiteralPath ''%s'' -Algorithm SHA256).Hash"'],literalPath);
    [exitCode,output]=system(command); digest=lower(strtrim(output));
    assert(exitCode==0&&~isempty(regexp(digest,'^[a-f0-9]{64}$','once')), ...
        'SHA256 failed: %s',output);
end
end

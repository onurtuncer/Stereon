param(
    [Parameter(Mandatory)] [string] $Dir,          # testcase directory; main.cpp is the last TU
    [string] $Compiler = "g162",
    [string] $Flags = "-std=c++26 -fmodules",
    [string[]] $Order = @(),                        # extra files, compiled before main.cpp
    [string[]] $Libraries = @(),                    # CE library ids as "id:version", e.g. stdexec:trunk
    [string] $Main = "main.cpp",                    # file sent as the main source
    [switch] $Execute,
    [switch] $CompileOnly,                          # -c instead of linking
    [string[]] $Tools = @(),                        # CE tools as "id:args", run on the output, e.g. "readelf:-sW"
    [switch] $Asm,                                  # return the (unfiltered) assembly of $Main
    [string] $AsmOut = ""                           # file to save the assembly to
)
$files = @()
foreach ($f in Get-ChildItem $Dir -File) {
    if ($f.Name -eq $Main -or $f.Name -like "*.txt" -or $f.Name -like "*.sh" -or $f.Name -like "*.md" -or $f.Name -like "*.s") { continue }
    $files += @{ filename = $f.Name; contents = [IO.File]::ReadAllText($f.FullName) }
}
$libs = @($Libraries | ForEach-Object { $p = $_.Split(":"); @{ id = $p[0]; version = $p[1] } })
$xs = ($Order | ForEach-Object { "-x c++ $_" }) -join " "
$userArgs = "$Flags $xs -x c++"
$filters = @{ binary = [bool](-not $CompileOnly -and -not $Asm); execute = [bool]$Execute }
if ($Asm) { $filters += @{ labels = $false; directives = $false; commentOnly = $true; libraryCode = $false; demangle = $false } }
$body = @{
    source = [IO.File]::ReadAllText((Join-Path $Dir $Main))
    options = @{
        userArguments = $userArgs
        compilerOptions = @{ executorRequest = [bool]$Execute; skipAsm = [bool](-not $Asm) }
        filters = $filters
        tools = @($Tools | ForEach-Object { $i = $_.IndexOf(":"); @{ id = $_.Substring(0, $i); args = $_.Substring($i + 1) } }); libraries = $libs
    }
    lang = "c++"
    files = $files
    allowStoreCodeDebug = $false
} | ConvertTo-Json -Depth 8
$r = Invoke-RestMethod -Method Post -Uri "https://godbolt.org/api/compiler/$Compiler/compile" -ContentType "application/json" -Headers @{Accept="application/json"} -Body $body
"== $Compiler  args: $userArgs  libs: $($Libraries -join ',')"
if ($Execute -and $r.buildResult) { $b = $r.buildResult; "build code: $($b.code)"; $b.stdout.text; $b.stderr.text; "run code: $($r.code)"; $r.stdout.text; $r.stderr.text }
else { "code: $($r.code)"; $r.stdout.text; $r.stderr.text }
foreach ($t in $r.tools) { "== tool $($t.id) code $($t.code)"; ($t.stdout | ForEach-Object { $_.text }) -join "`n"; ($t.stderr | ForEach-Object { $_.text }) -join "`n" }
if ($Asm) {
    $text = ($r.asm | ForEach-Object { $_.text }) -join "`n"
    if ($AsmOut) { [IO.File]::WriteAllText($AsmOut, $text); "asm lines: $($r.asm.Count) -> $AsmOut" } else { $text }
}

param(
    [Parameter(Mandatory)] [string] $Dir,          # testcase directory; main.cpp is the last TU
    [string] $Compiler = "g162",
    [string] $Flags = "-std=c++26 -fmodules",
    [string[]] $Order = @(),                        # extra files, compiled before main.cpp
    [switch] $Execute,
    [switch] $CompileOnly                           # no link; single-file testcases only (CE adds -S -o)
)
$files = @()
foreach ($f in Get-ChildItem $Dir -File) {
    if ($f.Name -eq "main.cpp" -or $f.Name -like "*.txt" -or $f.Name -like "*.sh" -or $f.Name -like "*.md") { continue }
    $files += @{ filename = $f.Name; contents = [IO.File]::ReadAllText($f.FullName) }
}
$xs = ($Order | ForEach-Object { "-x c++ $_" }) -join " "
$userArgs = "$Flags $xs -x c++"
$body = @{
    source = [IO.File]::ReadAllText((Join-Path $Dir "main.cpp"))
    options = @{
        userArguments = $userArgs
        compilerOptions = @{ executorRequest = [bool]$Execute; skipAsm = $true }
        filters = @{ binary = [bool](-not $CompileOnly); execute = [bool]$Execute }
        tools = @(); libraries = @()
    }
    lang = "c++"
    files = $files
    allowStoreCodeDebug = $false
} | ConvertTo-Json -Depth 8
$r = Invoke-RestMethod -Method Post -Uri "https://godbolt.org/api/compiler/$Compiler/compile" -ContentType "application/json" -Headers @{Accept="application/json"} -Body $body
"== $Compiler  args: $userArgs"
if ($Execute -and $r.buildResult) { $b = $r.buildResult; "build code: $($b.code)"; $b.stdout.text; $b.stderr.text; "run code: $($r.code)"; $r.stdout.text; $r.stderr.text }
else { "code: $($r.code)"; $r.stdout.text; $r.stderr.text }

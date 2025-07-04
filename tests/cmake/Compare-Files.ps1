<#
Compare-Files.ps1
Compares two text files while ignoring:
  - blank lines
  - differences in the amount of whitespace
Exit code: 0 = identical / 1 = different
#>

param(
  [Parameter(Mandatory)][string]$TestFile,
  [Parameter(Mandatory)][string]$OrigFile
)

function Normalize([string]$path) {
  Get-Content $path -Raw -Delimiter "`n" |
    ForEach-Object {
      $_ -replace "\r", ""      # CRLF → LF
    } | ForEach-Object {
      ($_ -replace "\s+", " ").Trim() # collapse whitespace & trim
    } | Where-Object { $_ -ne "" }    # drop blank lines
}

$test = Normalize $TestFile
$orig = Normalize $OrigFile

if ($test -ne $orig) {
  Write-Host "Difference detected between $TestFile and $OrigFile"
  Compare-Object $test $orig | Out-String | Write-Host
  exit 1
}
exit 0

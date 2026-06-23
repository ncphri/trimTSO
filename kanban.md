# trimTSO マルチスレッド化 Kanban

## タスク一覧

| 状態 | ID | 内容 |
|------|----|------|
| ✅ DONE | T1 | ストリーム入力 Reader 実装（gz/plain対応） |
| ✅ DONE | T2 | FastqProcessor チャンク処理＋スレッドプール |
| ✅ DONE | T3 | main() に -t オプション追加＋ループ構造 |
| ✅ DONE | T4 | 出力バッファリング削除（チャンク単位書き込み） |
| ✅ DONE | T5 | コンパイル＆動作確認 |
| ✅ DONE | T6 | 統計累積バグ修正（スレッドごとにゼロ初期化＋加算） |

## 修正内容 (T6)

**バグ:** 複数チャンクで累積 trimCount がスレッド数倍されていた
**原因:** 各スレッドの AdapterTrimmer が前チャンクまでの累積 trimCount を含んだコピーで初期化され、マージ時にスレッド間で重複加算されていた

**修正:**
1. `threadWorker` 内で `cleanAdapter` を作成（trimCount=0のコピー）→ スレッドは自分のトリムだけ計上
2. `processChunk` のマージロジックを `+=` に変更（累積に加算）

## 検証結果

| テスト | t=1 | t=2 | t=4 |
|--------|-----|-----|-----|
| Total trimmed reads | 2 | 2 | 2 ✅ |

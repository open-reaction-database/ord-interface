const {defineConfig} = require('@vue/cli-service')
const path = require('path')

module.exports = defineConfig({
    transpileDependencies: true,
    devServer: {
        proxy: {
            "^/api": {
                target: process.env.NODE_ENV === 'production' ? "https://open-reaction-database.org/client/" : "http://0.0.0.0:5000/client/",
                changeOrigin: true
            },
            "^/editor-api": {
                target: process.env.NODE_ENV === 'production' ? "https://open-reaction-database.org/editor/" : "http://0.0.0.0:5000/editor/",
                changeOrigin: true
            }
        }
    },
    chainWebpack: config => {
        config.plugin('copy')
            .tap(entries => {
                entries[0].patterns.push({
                    from: path.resolve(__dirname, 'src/ketcher/templates'),
                    to: path.resolve(__dirname, 'dist/templates'),
                    toType: 'dir',
                    noErrorOnMissing: false,
                    globOptions: {ignore: ['.DS_Store']},
                })
                return entries
            })
    }
})
